/* Render large, high-resolution Altadena maps from MODIS burned area and leaf area.
 *
 *   make maps
 *
 * Reads Cloud Optimized GeoTIFFs from the Planetary Computer STAC API and
 * composites them on an Esri World Topographic basemap. Output:
 *   output/altadena_mcd64a1_burn.png
 *   output/altadena_lai.png
 *   output/altadena_lai_summary.json
 */

#include <cairo.h>
#include <cairo-ft.h>
#include <curl/curl.h>
#include <ft2build.h>
#include FT_FREETYPE_H
#include <geotiff.h>
#include <geo_normalize.h>
#include <math.h>
#include <proj.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/stat.h>
#include <tiffio.h>
#include <xtiffio.h>

#define PANEL_W 3000
#define MARGIN 56
#define OVERLAY_ALPHA 0.78

static const double WEST = -118.20;
static const double SOUTH = 34.155;
static const double EAST = -118.03;
static const double NORTH = 34.265;
static const double MERC_R = 6378137.0;

static const char *STAC = "https://planetarycomputer.microsoft.com/api/stac/v1";

typedef struct {
    unsigned char *data;
    size_t len;
    size_t cap;
} Buf;

typedef struct {
    double r, g, b;
} RGB;

typedef struct {
    double x, y, w, h;
} Rect;

typedef struct {
    int w, h;
    double origin_x;
    double origin_y;
    double pixel_w;
    double pixel_h;
    float *v;
    PJ *to_merc;
} Grid;

typedef struct {
    const char *name;
    double lon, lat;
    int anchor;
} Place;

enum { ANCHOR_NE, ANCHOR_SE, ANCHOR_NW, ANCHOR_SW, ANCHOR_N, ANCHOR_S };

static const Place PLACES[] = {
    {"JPL", -118.171, 34.201, ANCHOR_NE},
    {"Farnsworth Park", -118.128, 34.203, ANCHOR_N},
    {"Fair Oaks / Loma Alta", -118.151, 34.193, ANCHOR_SW},
    {"Christmas Tree Lane", -118.131, 34.189, ANCHOR_S},
    {"Eaton Canyon", -118.096, 34.190, ANCHOR_NE},
    {"Kinneloa Mesa", -118.078, 34.185, ANCHOR_SE},
};

static const RGB YLGN[] = {
    {1.0, 1.0, 0.89803921568627454},
    {0.96862745098039216, 0.9882352941176471, 0.72549019607843135},
    {0.85098039215686272, 0.94117647058823528, 0.63921568627450975},
    {0.67843137254901964, 0.8666666666666667, 0.55686274509803924},
    {0.47058823529411764, 0.77647058823529413, 0.47450980392156861},
    {0.25490196078431371, 0.6705882352941176, 0.36470588235294116},
    {0.13725490196078433, 0.51764705882352946, 0.2627450980392157},
    {0.0, 0.40784313725490196, 0.21568627450980393},
    {0.0, 0.27058823529411763, 0.16078431372549021},
};

static const RGB RDYLGN[] = {
    {0.6470588235294118, 0.0, 0.14901960784313725},
    {0.84313725490196079, 0.18823529411764706, 0.15294117647058825},
    {0.95686274509803926, 0.42745098039215684, 0.2627450980392157},
    {0.99215686274509807, 0.68235294117647061, 0.38039215686274508},
    {0.99607843137254903, 0.8784313725490196, 0.54509803921568623},
    {1.0, 1.0, 0.74901960784313726},
    {0.85098039215686272, 0.93725490196078431, 0.54509803921568623},
    {0.65098039215686276, 0.85098039215686272, 0.41568627450980394},
    {0.4, 0.74117647058823533, 0.38823529411764707},
    {0.10196078431372549, 0.59607843137254901, 0.31372549019607843},
    {0.0, 0.40784313725490196, 0.21568627450980393},
};

static const RGB BURN_COLORS[] = {
    {254 / 255.0, 224 / 255.0, 139 / 255.0},
    {253 / 255.0, 174 / 255.0, 97 / 255.0},
    {244 / 255.0, 109 / 255.0, 67 / 255.0},
    {215 / 255.0, 48 / 255.0, 39 / 255.0},
    {165 / 255.0, 0 / 255.0, 38 / 255.0},
    {103 / 255.0, 0 / 255.0, 31 / 255.0},
};

static PJ_CONTEXT *g_proj;
static char *g_token;

static void die(const char *msg) {
    fprintf(stderr, "error: %s\n", msg);
    exit(1);
}

static void *xmalloc(size_t n) {
    void *p = malloc(n);
    if (!p) die("out of memory");
    return p;
}

static char *xstrdup(const char *s) {
    char *p = xmalloc(strlen(s) + 1);
    memcpy(p, s, strlen(s) + 1);
    return p;
}

static int finite_f(float v) { return v == v; }

static void lonlat_merc(double lon, double lat, double *x, double *y) {
    double phi = lat * M_PI / 180.0;
    *x = MERC_R * lon * M_PI / 180.0;
    *y = MERC_R * log(tan(M_PI * 0.25 + phi * 0.5));
}

static RGB lerp_stops(const RGB *stops, int n, double t) {
    if (t < 0) t = 0;
    if (t > 1) t = 1;
    double x = t * (n - 1);
    int i = (int)x;
    if (i >= n - 1) return stops[n - 1];
    double f = x - i;
    RGB c = {
        stops[i].r + (stops[i + 1].r - stops[i].r) * f,
        stops[i].g + (stops[i + 1].g - stops[i].g) * f,
        stops[i].b + (stops[i + 1].b - stops[i].b) * f,
    };
    return c;
}

static int burn_bin(float day) {
    const double edges[] = {6, 8, 9, 10, 11, 12, 14};
    if (day < edges[0]) return 0;
    for (int i = 0; i < 6; i++) {
        if (day >= edges[i] && day < edges[i + 1]) return i;
    }
    return 5;
}

static size_t write_mem(char *ptr, size_t size, size_t nmemb, void *ud) {
    size_t n = size * nmemb;
    Buf *b = ud;
    if (b->len + n + 1 > b->cap) {
        size_t cap = b->cap ? b->cap : 65536;
        while (cap < b->len + n + 1) cap *= 2;
        unsigned char *p = realloc(b->data, cap);
        if (!p) return 0;
        b->data = p;
        b->cap = cap;
    }
    memcpy(b->data + b->len, ptr, n);
    b->len += n;
    b->data[b->len] = 0;
    return n;
}

static Buf http_get(CURL *curl, const char *url, long *code) {
    Buf b = {0};
    char err[CURL_ERROR_SIZE];
    err[0] = 0;
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, write_mem);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, &b);
    curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, err);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 180L);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "altadena-maps/1.0");
    CURLcode rc = curl_easy_perform(curl);
    if (rc != CURLE_OK) {
        fprintf(stderr, "curl: %s\n", err[0] ? err : curl_easy_strerror(rc));
        free(b.data);
        b.data = NULL;
        b.len = 0;
        if (code) *code = 0;
        return b;
    }
    long status = 0;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &status);
    if (code) *code = status;
    if (status != 200) {
        free(b.data);
        b.data = NULL;
        b.len = 0;
    }
    return b;
}

static int http_download(CURL *curl, const char *url, const char *path) {
    FILE *f = fopen(path, "wb");
    if (!f) die("could not write cache file");
    char err[CURL_ERROR_SIZE];
    err[0] = 0;
    curl_easy_setopt(curl, CURLOPT_URL, url);
    curl_easy_setopt(curl, CURLOPT_WRITEFUNCTION, NULL);
    curl_easy_setopt(curl, CURLOPT_WRITEDATA, f);
    curl_easy_setopt(curl, CURLOPT_ERRORBUFFER, err);
    curl_easy_setopt(curl, CURLOPT_FOLLOWLOCATION, 1L);
    curl_easy_setopt(curl, CURLOPT_TIMEOUT, 180L);
    curl_easy_setopt(curl, CURLOPT_USERAGENT, "altadena-maps/1.0");
    CURLcode rc = curl_easy_perform(curl);
    fclose(f);
    if (rc != CURLE_OK) {
        fprintf(stderr, "curl: %s\n", err[0] ? err : curl_easy_strerror(rc));
        remove(path);
        return 0;
    }
    long status = 0;
    long redirects = 0;
    char *effective = NULL;
    curl_easy_getinfo(curl, CURLINFO_RESPONSE_CODE, &status);
    curl_easy_getinfo(curl, CURLINFO_REDIRECT_COUNT, &redirects);
    curl_easy_getinfo(curl, CURLINFO_EFFECTIVE_URL, &effective);
    if (status != 200) {
        fprintf(stderr, "download HTTP %ld (redirects %ld, query %s)\n", status, redirects,
                effective && strstr(effective, "sig=") ? "kept" : "missing");
        FILE *errf = fopen(path, "rb");
        if (errf) {
            char snippet[240];
            size_t n = fread(snippet, 1, sizeof snippet - 1, errf);
            snippet[n] = 0;
            fprintf(stderr, "response: %s\n", snippet);
            fclose(errf);
        }
        remove(path);
        return 0;
    }
    return 1;
}

static int hex_nibble(char h) {
    if (h >= '0' && h <= '9') return h - '0';
    if (h >= 'a' && h <= 'f') return h - 'a' + 10;
    if (h >= 'A' && h <= 'F') return h - 'A' + 10;
    return -1;
}

/* Copy a JSON string, turning \\u0026 and similar escapes back into characters. */
static char *json_unescape(const char *start, const char *end) {
    char *out = xmalloc((size_t)(end - start) + 1);
    size_t j = 0;
    for (const char *p = start; p < end; p++) {
        if (*p != '\\' || p + 1 >= end) {
            out[j++] = *p;
            continue;
        }
        char e = *++p;
        if (e == 'u' && p + 4 < end) {
            int cp = 0;
            int ok = 1;
            for (int k = 0; k < 4; k++) {
                int n = hex_nibble(*++p);
                if (n < 0) ok = 0;
                cp = (cp << 4) | n;
            }
            out[j++] = ok && cp >= 0 && cp < 128 ? (char)cp : '?';
        } else if (e == 'n') {
            out[j++] = '\n';
        } else if (e == 't') {
            out[j++] = '\t';
        } else if (e == 'r') {
            out[j++] = '\r';
        } else {
            out[j++] = e;
        }
    }
    out[j] = 0;
    return out;
}

static char *json_quoted(const char *json, const char *key) {
    char pat[160];
    snprintf(pat, sizeof pat, "\"%s\"", key);
    const char *p = strstr(json, pat);
    if (!p) return NULL;
    p += strlen(pat);
    while (*p && *p != '"') p++;
    if (*p != '"') return NULL;
    p++;
    const char *end = p;
    while (*end && *end != '"') {
        if (*end == '\\' && end[1]) end += 2;
        else end++;
    }
    if (*end != '"') return NULL;
    return json_unescape(p, end);
}

static char *sas_token(CURL *curl) {
    if (g_token) return g_token;
    long code = 0;
    Buf body = http_get(
        curl,
        "https://planetarycomputer.microsoft.com/api/sas/v1/token/modiseuwest/modis-061-cogs",
        &code);
    if (!body.data) {
        fprintf(stderr, "SAS HTTP %ld\n", code);
        die("could not get a Planetary Computer SAS token");
    }
    g_token = json_quoted((char *)body.data, "token");
    free(body.data);
    if (!g_token) die("SAS token response had no token");
    return g_token;
}

static char *asset_href(CURL *curl, const char *collection, const char *item, const char *asset) {
    char url[768];
    snprintf(url, sizeof url, "%s/collections/%s/items/%s", STAC, collection, item);
    long code = 0;
    Buf body = http_get(curl, url, &code);
    if (!body.data) die("could not fetch STAC item");
    char pat[160];
    snprintf(pat, sizeof pat, "\"%s\"", asset);
    const char *p = strstr((char *)body.data, pat);
    if (!p) die("asset not found on STAC item");
    char *href = json_quoted(p, "href");
    free(body.data);
    if (!href) die("asset href missing");
    return href;
}

static int file_ok(const char *path) {
    struct stat st;
    return stat(path, &st) == 0 && st.st_size > 1000;
}

static uint16_t *read_band(TIFF *tif, uint32_t w, uint32_t h, uint16_t bps) {
    if (bps != 8 && bps != 16) die("expected an 8-bit or 16-bit MODIS raster");
    uint16_t *out = xmalloc((size_t)w * h * sizeof(uint16_t));
    if (TIFFIsTiled(tif)) {
        uint32_t tw = 0, th = 0;
        TIFFGetField(tif, TIFFTAG_TILEWIDTH, &tw);
        TIFFGetField(tif, TIFFTAG_TILELENGTH, &th);
        tmsize_t tsize = TIFFTileSize(tif);
        unsigned char *tile = xmalloc((size_t)tsize);
        for (uint32_t row = 0; row < h; row += th) {
            for (uint32_t col = 0; col < w; col += tw) {
                if (TIFFReadTile(tif, tile, col, row, 0, 0) < 0) die("could not read a GeoTIFF tile");
                uint32_t rows = row + th > h ? h - row : th;
                uint32_t cols = col + tw > w ? w - col : tw;
                for (uint32_t r = 0; r < rows; r++) {
                    for (uint32_t c = 0; c < cols; c++) {
                        uint16_t v;
                        if (bps == 8) v = tile[(size_t)r * tw + c];
                        else v = ((uint16_t *)tile)[(size_t)r * tw + c];
                        out[(size_t)(row + r) * w + (col + c)] = v;
                    }
                }
            }
        }
        free(tile);
    } else {
        tmsize_t sl = TIFFScanlineSize(tif);
        unsigned char *buf = xmalloc((size_t)sl);
        for (uint32_t row = 0; row < h; row++) {
            if (TIFFReadScanline(tif, buf, row, 0) < 0) die("could not read a GeoTIFF scanline");
            for (uint32_t c = 0; c < w; c++) {
                uint16_t v = bps == 8 ? buf[c] : ((uint16_t *)buf)[c];
                out[(size_t)row * w + c] = v;
            }
        }
        free(buf);
    }
    return out;
}

static char *crs_from_proj4(const char *proj4) {
    if (strstr(proj4, "+type=crs")) return xstrdup(proj4);
    size_t n = strlen(proj4);
    char *s = xmalloc(n + 12);
    memcpy(s, proj4, n);
    memcpy(s + n, " +type=crs", 11);
    return s;
}

static PJ *crs_op(const char *src, const char *dst, int normalize) {
    PJ *op = proj_create_crs_to_crs(g_proj, src, dst, NULL);
    if (!op) die(proj_context_errno_string(g_proj, proj_context_errno(g_proj)));
    if (!normalize) return op;
    PJ *vis = proj_normalize_for_visualization(g_proj, op);
    proj_destroy(op);
    if (!vis) die("could not set lon/lat axis order");
    return vis;
}

static Grid *load_grid(CURL *curl, const char *collection, const char *item, const char *asset, int kind_lai) {
    char cache[512];
    snprintf(cache, sizeof cache, "cache/%s_%s.tif", item, asset);
    for (char *p = cache + 6; *p; p++) {
        if (*p == '/') *p = '_';
    }
    if (!file_ok(cache)) {
        char *href = asset_href(curl, collection, item, asset);
        const char *token = sas_token(curl);
        size_t href_n = strlen(href);
        size_t token_n = strlen(token);
        char *url = xmalloc(href_n + token_n + 2);
        memcpy(url, href, href_n);
        url[href_n] = '?';
        memcpy(url + href_n + 1, token, token_n);
        url[href_n + token_n + 1] = 0;
        fprintf(stderr, "downloading %s\n", href);
        if (!http_download(curl, url, cache)) die("GeoTIFF download failed");
        free(url);
        free(href);
    } else {
        fprintf(stderr, "using cached %s\n", cache);
    }

    TIFF *tif = XTIFFOpen(cache, "r");
    if (!tif) die("could not open GeoTIFF");
    GTIF *gt = GTIFNew(tif);
    if (!gt) die("could not read GeoTIFF keys");
    GTIFDefn defn;
    if (!GTIFGetDefn(gt, &defn)) die("could not decode GeoTIFF projection");
    char *proj4 = GTIFGetProj4Defn(&defn);
    if (!proj4) die("could not build a PROJ string for this raster");

    double x0 = 0, y0 = 0, x1 = 1, y1 = 0, x2 = 0, y2 = 1;
    if (!GTIFImageToPCS(gt, &x0, &y0) || !GTIFImageToPCS(gt, &x1, &y1) || !GTIFImageToPCS(gt, &x2, &y2))
        die("could not map pixel corners");
    double pixel_w = x1 - x0;
    double pixel_h = y0 - y2;
    if (pixel_w <= 0 || pixel_h <= 0) die("unexpected pixel size");

    uint32_t iw = 0, ih = 0;
    uint16_t bps = 8;
    TIFFGetField(tif, TIFFTAG_IMAGEWIDTH, &iw);
    TIFFGetField(tif, TIFFTAG_IMAGELENGTH, &ih);
    TIFFGetField(tif, TIFFTAG_BITSPERSAMPLE, &bps);
    uint16_t *raw = read_band(tif, iw, ih, bps);
    GTIFFree(gt);
    XTIFFClose(tif);

    char *src = crs_from_proj4(proj4);
    free(proj4);
    PJ *from_ll = crs_op("EPSG:4326", src, 1);
    PJ *to_merc = crs_op(src, "EPSG:3857", 0);
    free(src);

    PJ_COORD sw = proj_trans(from_ll, PJ_FWD, proj_coord(WEST, SOUTH, 0, 0));
    PJ_COORD ne = proj_trans(from_ll, PJ_FWD, proj_coord(EAST, NORTH, 0, 0));
    if (sw.xy.x == HUGE_VAL || ne.xy.x == HUGE_VAL) die("could not project the map bounds");
    double left = fmin(sw.xy.x, ne.xy.x);
    double right = fmax(sw.xy.x, ne.xy.x);
    double bottom = fmin(sw.xy.y, ne.xy.y);
    double top = fmax(sw.xy.y, ne.xy.y);
    double col_off = (left - x0) / pixel_w;
    double row_off = (y0 - top) / pixel_h;
    double cols = (right - left) / pixel_w;
    double rows = (top - bottom) / pixel_h;
    int c0 = (int)lround(col_off);
    int r0 = (int)lround(row_off);
    int c1 = (int)lround(col_off + cols);
    int r1 = (int)lround(row_off + rows);
    if (c0 < 0) c0 = 0;
    if (r0 < 0) r0 = 0;
    if (c1 > (int)iw) c1 = (int)iw;
    if (r1 > (int)ih) r1 = (int)ih;
    if (c1 <= c0 || r1 <= r0) die("map window missed the MODIS tile");

    Grid *g = calloc(1, sizeof(Grid));
    if (!g) die("out of memory");
    g->w = c1 - c0;
    g->h = r1 - r0;
    g->origin_x = x0 + c0 * pixel_w;
    g->origin_y = y0 - r0 * pixel_h;
    g->pixel_w = pixel_w;
    g->pixel_h = pixel_h;
    g->to_merc = to_merc;
    g->v = xmalloc((size_t)g->w * g->h * sizeof(float));
    for (int r = 0; r < g->h; r++) {
        for (int c = 0; c < g->w; c++) {
            uint16_t sample = raw[(size_t)(r0 + r) * iw + (c0 + c)];
            float v = NAN;
            if (kind_lai) {
                if (sample <= 100) v = sample * 0.1f;
            } else if (sample >= 1 && sample <= 366) {
                v = sample;
            }
            g->v[(size_t)r * g->w + c] = v;
        }
    }
    free(raw);
    proj_destroy(from_ll);
    fprintf(stderr, "  window %d x %d px, pixel %.1f m\n", g->w, g->h, pixel_w);
    if (g->w > 500 || g->h > 500 || g->w < 8 || g->h < 8)
        die("map window size is implausible; the projection did not land on Altadena");
    return g;
}

static void free_grid(Grid *g) {
    if (!g) return;
    proj_destroy(g->to_merc);
    free(g->v);
    free(g);
}

static void cell_center(const Grid *g, int c, int r, double *x, double *y) {
    *x = g->origin_x + (c + 0.5) * g->pixel_w;
    *y = g->origin_y - (r + 0.5) * g->pixel_h;
}

static float sample_nn(const Grid *g, double x, double y) {
    int c = (int)floor((x - g->origin_x) / g->pixel_w);
    int r = (int)floor((g->origin_y - y) / g->pixel_h);
    if (c < 0 || r < 0 || c >= g->w || r >= g->h) return NAN;
    return g->v[(size_t)r * g->w + c];
}

static void project_corner(const Grid *g, int c, int r, double *mx, double *my) {
    double x = g->origin_x + c * g->pixel_w;
    double y = g->origin_y - r * g->pixel_h;
    PJ_COORD o = proj_trans(g->to_merc, PJ_FWD, proj_coord(x, y, 0, 0));
    if (o.xy.x == HUGE_VAL) die("could not project a MODIS pixel");
    *mx = o.xy.x;
    *my = o.xy.y;
}

typedef struct {
    double left, right, bottom, top;
    Rect panel;
} View;

static void merc_to_screen(const View *v, double mx, double my, double *sx, double *sy) {
    *sx = v->panel.x + (mx - v->left) / (v->right - v->left) * v->panel.w;
    *sy = v->panel.y + (v->top - my) / (v->top - v->bottom) * v->panel.h;
}

typedef struct {
    double x0, y0, x1, y1;
} Box;

static Box *screen_boxes(const Grid *g, const View *v) {
    int cw = g->w + 1, ch = g->h + 1;
    double *sx = xmalloc((size_t)cw * ch * sizeof(double));
    double *sy = xmalloc((size_t)cw * ch * sizeof(double));
    for (int r = 0; r < ch; r++) {
        for (int c = 0; c < cw; c++) {
            double mx, my;
            project_corner(g, c, r, &mx, &my);
            merc_to_screen(v, mx, my, &sx[(size_t)r * cw + c], &sy[(size_t)r * cw + c]);
        }
    }
    Box *boxes = xmalloc((size_t)g->w * g->h * sizeof(Box));
    for (int r = 0; r < g->h; r++) {
        for (int c = 0; c < g->w; c++) {
            double x[4] = {
                sx[(size_t)r * cw + c],
                sx[(size_t)r * cw + c + 1],
                sx[(size_t)(r + 1) * cw + c],
                sx[(size_t)(r + 1) * cw + c + 1],
            };
            double y[4] = {
                sy[(size_t)r * cw + c],
                sy[(size_t)r * cw + c + 1],
                sy[(size_t)(r + 1) * cw + c],
                sy[(size_t)(r + 1) * cw + c + 1],
            };
            double x0 = x[0], x1 = x[0], y0 = y[0], y1 = y[0];
            for (int k = 1; k < 4; k++) {
                if (x[k] < x0) x0 = x[k];
                if (x[k] > x1) x1 = x[k];
                if (y[k] < y0) y0 = y[k];
                if (y[k] > y1) y1 = y[k];
            }
            boxes[(size_t)r * g->w + c] = (Box){x0 - 0.8, y0 - 0.8, x1 + 0.8, y1 + 0.8};
        }
    }
    free(sx);
    free(sy);
    return boxes;
}

static void fill_cells(cairo_t *cr, const Grid *g, const Box *boxes, RGB (*color_of)(float, void *), void *ud) {
    for (int r = 0; r < g->h; r++) {
        for (int c = 0; c < g->w; c++) {
            float v = g->v[(size_t)r * g->w + c];
            if (!finite_f(v)) continue;
            RGB col = color_of(v, ud);
            Box b = boxes[(size_t)r * g->w + c];
            cairo_set_source_rgba(cr, col.r, col.g, col.b, OVERLAY_ALPHA);
            cairo_rectangle(cr, b.x0, b.y0, b.x1 - b.x0, b.y1 - b.y0);
            cairo_fill(cr);
        }
    }
}

static RGB color_burn(float day, void *ud) {
    (void)ud;
    return BURN_COLORS[burn_bin(day)];
}

static RGB color_lai(float lai, void *ud) {
    (void)ud;
    double t = lai / 2.5;
    return lerp_stops(YLGN, 9, t);
}

static RGB color_change(float d, void *ud) {
    (void)ud;
    double t = (d + 1.5) / 3.0;
    return lerp_stops(RDYLGN, 11, t);
}

static cairo_status_t read_png_buf(void *closure, unsigned char *data, unsigned int length) {
    Buf *b = closure;
    if (b->cap + length > b->len) return CAIRO_STATUS_READ_ERROR;
    memcpy(data, b->data + b->cap, length);
    b->cap += length;
    return CAIRO_STATUS_SUCCESS;
}

static cairo_surface_t *fetch_basemap(CURL *curl, double left, double bottom, double right, double top, int w, int h) {
    char url[640];
    snprintf(url, sizeof url,
             "https://server.arcgisonline.com/ArcGIS/rest/services/World_Topo_Map/MapServer/export"
             "?bbox=%.3f,%.3f,%.3f,%.3f&bboxSR=3857&imageSR=3857&size=%d,%d&format=png32&f=image",
             left, bottom, right, top, w, h);
    fprintf(stderr, "fetching basemap %d x %d\n", w, h);
    long code = 0;
    Buf body = {0};
    for (int attempt = 0; attempt < 3 && !body.data; attempt++) {
        body = http_get(curl, url, &code);
        if (!body.data) fprintf(stderr, "basemap attempt %d failed (HTTP %ld)\n", attempt + 1, code);
    }
    if (!body.data) die("could not download the basemap");
    if (body.len < 8 || memcmp(body.data, "\x89PNG", 4) != 0) die("basemap response was not a PNG");
    Buf reader = body;
    reader.cap = 0;
    cairo_surface_t *img = cairo_image_surface_create_from_png_stream(read_png_buf, &reader);
    free(body.data);
    if (cairo_surface_status(img) != CAIRO_STATUS_SUCCESS) die("could not decode the basemap");
    return img;
}

static void paint_basemap(cairo_t *cr, cairo_surface_t *base, Rect panel) {
    cairo_save(cr);
    cairo_rectangle(cr, panel.x, panel.y, panel.w, panel.h);
    cairo_clip(cr);
    cairo_translate(cr, panel.x, panel.y);
    double bw = cairo_image_surface_get_width(base);
    double bh = cairo_image_surface_get_height(base);
    cairo_scale(cr, panel.w / bw, panel.h / bh);
    cairo_set_source_surface(cr, base, 0, 0);
    cairo_pattern_set_filter(cairo_get_source(cr), CAIRO_FILTER_BILINEAR);
    cairo_paint(cr);
    cairo_restore(cr);
}

static void round_rect(cairo_t *cr, double x, double y, double w, double h, double r) {
    if (r > w / 2) r = w / 2;
    if (r > h / 2) r = h / 2;
    cairo_new_sub_path(cr);
    cairo_arc(cr, x + w - r, y + r, r, -M_PI / 2, 0);
    cairo_arc(cr, x + w - r, y + h - r, r, 0, M_PI / 2);
    cairo_arc(cr, x + r, y + h - r, r, M_PI / 2, M_PI);
    cairo_arc(cr, x + r, y + r, r, M_PI, 3 * M_PI / 2);
    cairo_close_path(cr);
}

static void draw_label(cairo_t *cr, cairo_font_face_t *face, double size, double px, double py,
                       const char *text, int anchor) {
    cairo_set_font_face(cr, face);
    cairo_set_font_size(cr, size);
    cairo_text_extents_t e;
    cairo_text_extents(cr, text, &e);
    double pad = size * 0.22;
    double tw = e.width + pad * 2;
    double th = size + pad * 2;
    double gap = size * 0.28;
    double x = px, y = py;
    if (anchor == ANCHOR_NE || anchor == ANCHOR_SE) x = px + gap;
    if (anchor == ANCHOR_NW || anchor == ANCHOR_SW) x = px - gap - tw;
    if (anchor == ANCHOR_N || anchor == ANCHOR_S) x = px - tw / 2;
    if (anchor == ANCHOR_NE || anchor == ANCHOR_NW || anchor == ANCHOR_N) y = py - gap - th;
    if (anchor == ANCHOR_SE || anchor == ANCHOR_SW || anchor == ANCHOR_S) y = py + gap;
    cairo_set_source_rgba(cr, 1, 1, 1, 0.84);
    round_rect(cr, x, y, tw, th, size * 0.18);
    cairo_fill(cr);
    cairo_set_source_rgb(cr, 0.10, 0.10, 0.10);
    cairo_move_to(cr, x + pad - e.x_bearing, y + (th - e.height) / 2 - e.y_bearing);
    cairo_show_text(cr, text);
}

static void draw_places(cairo_t *cr, cairo_font_face_t *face, const View *v) {
    double size = v->panel.w * 0.0115;
    for (size_t i = 0; i < sizeof PLACES / sizeof PLACES[0]; i++) {
        double mx, my, sx, sy;
        lonlat_merc(PLACES[i].lon, PLACES[i].lat, &mx, &my);
        merc_to_screen(v, mx, my, &sx, &sy);
        cairo_set_source_rgba(cr, 1, 1, 1, 0.95);
        cairo_arc(cr, sx, sy, size * 0.22, 0, 2 * M_PI);
        cairo_fill(cr);
        cairo_set_source_rgb(cr, 0.08, 0.08, 0.08);
        cairo_arc(cr, sx, sy, size * 0.13, 0, 2 * M_PI);
        cairo_fill(cr);
        draw_label(cr, face, size, sx, sy, PLACES[i].name, PLACES[i].anchor);
    }
}

static void draw_scale(cairo_t *cr, cairo_font_face_t *face, const View *v) {
    double lat = (SOUTH + NORTH) * 0.5 * M_PI / 180.0;
    double merc_for_2km = 2000.0 / cos(lat);
    double bar = merc_for_2km / (v->right - v->left) * v->panel.w;
    double x = v->panel.x + 28;
    double y = v->panel.y + v->panel.h - 36;
    cairo_set_font_face(cr, face);
    cairo_set_font_size(cr, 22);
    cairo_text_extents_t e;
    cairo_text_extents(cr, "2 km", &e);
    double baseline = y - 14;
    double text_top = baseline + e.y_bearing;
    double bw = bar + 28;
    if (e.width + 20 > bw) bw = e.width + 20;
    cairo_set_source_rgba(cr, 1, 1, 1, 0.9);
    round_rect(cr, x - 12, text_top - 8, bw, (y + 16) - (text_top - 8), 6);
    cairo_fill(cr);
    cairo_set_source_rgb(cr, 0.12, 0.12, 0.12);
    cairo_set_line_width(cr, 3);
    cairo_move_to(cr, x, y);
    cairo_line_to(cr, x + bar, y);
    cairo_move_to(cr, x, y - 7);
    cairo_line_to(cr, x, y + 7);
    cairo_move_to(cr, x + bar, y - 7);
    cairo_line_to(cr, x + bar, y + 7);
    cairo_stroke(cr);
    cairo_move_to(cr, x, baseline);
    cairo_show_text(cr, "2 km");
}

static void draw_centered(cairo_t *cr, cairo_font_face_t *face, double size, double cx, double baseline,
                          const char *text, RGB color) {
    cairo_set_font_face(cr, face);
    cairo_set_font_size(cr, size);
    cairo_set_source_rgb(cr, color.r, color.g, color.b);
    cairo_text_extents_t e;
    cairo_text_extents(cr, text, &e);
    cairo_move_to(cr, cx - (e.width / 2 + e.x_bearing), baseline);
    cairo_show_text(cr, text);
}

static void stroke_burn_edge(cairo_t *cr, const Grid *g, const View *v, const unsigned char *mask) {
    cairo_new_path(cr);
    for (int r = 0; r < g->h; r++) {
        for (int c = 0; c < g->w; c++) {
            if (!mask[(size_t)r * g->w + c]) continue;
            int left = c == 0 || !mask[(size_t)r * g->w + (c - 1)];
            int right = c + 1 == g->w || !mask[(size_t)r * g->w + (c + 1)];
            int up = r == 0 || !mask[(size_t)(r - 1) * g->w + c];
            int down = r + 1 == g->h || !mask[(size_t)(r + 1) * g->w + c];
            double mx[2], my[2], sx, sy;
            if (up) {
                project_corner(g, c, r, &mx[0], &my[0]);
                project_corner(g, c + 1, r, &mx[1], &my[1]);
                merc_to_screen(v, mx[0], my[0], &sx, &sy);
                cairo_move_to(cr, sx, sy);
                merc_to_screen(v, mx[1], my[1], &sx, &sy);
                cairo_line_to(cr, sx, sy);
            }
            if (down) {
                project_corner(g, c, r + 1, &mx[0], &my[0]);
                project_corner(g, c + 1, r + 1, &mx[1], &my[1]);
                merc_to_screen(v, mx[0], my[0], &sx, &sy);
                cairo_move_to(cr, sx, sy);
                merc_to_screen(v, mx[1], my[1], &sx, &sy);
                cairo_line_to(cr, sx, sy);
            }
            if (left) {
                project_corner(g, c, r, &mx[0], &my[0]);
                project_corner(g, c, r + 1, &mx[1], &my[1]);
                merc_to_screen(v, mx[0], my[0], &sx, &sy);
                cairo_move_to(cr, sx, sy);
                merc_to_screen(v, mx[1], my[1], &sx, &sy);
                cairo_line_to(cr, sx, sy);
            }
            if (right) {
                project_corner(g, c + 1, r, &mx[0], &my[0]);
                project_corner(g, c + 1, r + 1, &mx[1], &my[1]);
                merc_to_screen(v, mx[0], my[0], &sx, &sy);
                cairo_move_to(cr, sx, sy);
                merc_to_screen(v, mx[1], my[1], &sx, &sy);
                cairo_line_to(cr, sx, sy);
            }
        }
    }
}

static void draw_outline(cairo_t *cr, const Grid *g, const View *v, const unsigned char *mask) {
    cairo_save(cr);
    cairo_rectangle(cr, v->panel.x, v->panel.y, v->panel.w, v->panel.h);
    cairo_clip(cr);
    cairo_set_line_cap(cr, CAIRO_LINE_CAP_ROUND);
    cairo_set_line_join(cr, CAIRO_LINE_JOIN_ROUND);
    stroke_burn_edge(cr, g, v, mask);
    cairo_set_source_rgba(cr, 1, 1, 1, 0.95);
    cairo_set_line_width(cr, 8);
    cairo_stroke_preserve(cr);
    cairo_set_source_rgb(cr, 0.08, 0.08, 0.08);
    cairo_set_line_width(cr, 3);
    cairo_stroke(cr);
    cairo_restore(cr);
}

static void draw_panel_frame(cairo_t *cr, Rect panel) {
    cairo_set_source_rgba(cr, 0, 0, 0, 0.45);
    cairo_set_line_width(cr, 2);
    cairo_rectangle(cr, panel.x + 1, panel.y + 1, panel.w - 2, panel.h - 2);
    cairo_stroke(cr);
}

static void colorbar_continuous(cairo_t *cr, cairo_font_face_t *face, Rect bar, const RGB *stops, int nstops,
                                double vmin, double vmax, const char *title, const char *fmt) {
    cairo_set_font_face(cr, face);
    cairo_set_font_size(cr, 22);
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    cairo_text_extents_t e;
    cairo_text_extents(cr, title, &e);
    cairo_move_to(cr, bar.x - e.x_bearing, bar.y - 34);
    cairo_show_text(cr, title);

    int steps = (int)bar.h;
    for (int i = 0; i < steps; i++) {
        double t = 1.0 - (i + 0.5) / steps;
        RGB c = lerp_stops(stops, nstops, t);
        cairo_set_source_rgb(cr, c.r, c.g, c.b);
        cairo_rectangle(cr, bar.x, bar.y + i, bar.w, 1.2);
        cairo_fill(cr);
    }
    cairo_set_source_rgb(cr, 0.15, 0.15, 0.15);
    cairo_set_line_width(cr, 1.5);
    cairo_rectangle(cr, bar.x, bar.y, bar.w, bar.h);
    cairo_stroke(cr);

    cairo_set_font_size(cr, 20);
    int nt = 6;
    for (int i = 0; i < nt; i++) {
        double val = vmin + (vmax - vmin) * i / (nt - 1);
        double t = (val - vmin) / (vmax - vmin);
        double y = bar.y + bar.h - t * bar.h;
        cairo_set_source_rgb(cr, 0.15, 0.15, 0.15);
        cairo_set_line_width(cr, 1.5);
        cairo_move_to(cr, bar.x + bar.w, y);
        cairo_line_to(cr, bar.x + bar.w + 8, y);
        cairo_stroke(cr);
        char lab[32];
        snprintf(lab, sizeof lab, fmt, val);
        cairo_text_extents(cr, lab, &e);
        double base = y - (e.y_bearing + e.height / 2.0);
        double ink_top = base + e.y_bearing;
        double ink_bot = ink_top + e.height;
        if (ink_top < bar.y) base += bar.y - ink_top;
        if (ink_bot > bar.y + bar.h) base -= ink_bot - (bar.y + bar.h);
        cairo_move_to(cr, bar.x + bar.w + 12, base);
        cairo_show_text(cr, lab);
    }
}

static void colorbar_burn(cairo_t *cr, cairo_font_face_t *face, Rect bar) {
    cairo_set_font_face(cr, face);
    cairo_set_font_size(cr, 22);
    cairo_set_source_rgb(cr, 0.1, 0.1, 0.1);
    const char *title = "Burn day of 2025";
    cairo_text_extents_t e;
    cairo_text_extents(cr, title, &e);
    cairo_move_to(cr, bar.x - e.x_bearing, bar.y - 34);
    cairo_show_text(cr, title);
    const double edges[] = {6, 8, 9, 10, 11, 12, 14};
    double span = edges[6] - edges[0];
    for (int i = 0; i < 6; i++) {
        double y_top = bar.y + bar.h * (1.0 - (edges[i + 1] - edges[0]) / span);
        double y_bot = bar.y + bar.h * (1.0 - (edges[i] - edges[0]) / span);
        RGB c = BURN_COLORS[i];
        cairo_set_source_rgb(cr, c.r, c.g, c.b);
        cairo_rectangle(cr, bar.x, y_top, bar.w, y_bot - y_top + 0.5);
        cairo_fill(cr);
    }
    cairo_set_source_rgb(cr, 0.15, 0.15, 0.15);
    cairo_set_line_width(cr, 1.5);
    cairo_rectangle(cr, bar.x, bar.y, bar.w, bar.h);
    cairo_stroke(cr);
    cairo_set_font_size(cr, 20);
    for (int i = 0; i < 7; i++) {
        double t = (edges[i] - 6) / 8.0;
        double y = bar.y + bar.h - t * bar.h;
        cairo_move_to(cr, bar.x + bar.w, y);
        cairo_line_to(cr, bar.x + bar.w + 8, y);
        cairo_stroke(cr);
        char lab[16];
        snprintf(lab, sizeof lab, "%.0f", edges[i]);
        cairo_text_extents(cr, lab, &e);
        double base = y - (e.y_bearing + e.height / 2.0);
        double ink_top = base + e.y_bearing;
        double ink_bot = ink_top + e.height;
        if (ink_top < bar.y) base += bar.y - ink_top;
        if (ink_bot > bar.y + bar.h) base -= ink_bot - (bar.y + bar.h);
        cairo_move_to(cr, bar.x + bar.w + 12, base);
        cairo_show_text(cr, lab);
    }
}

static void draw_raster(cairo_t *cr, const Grid *g, const View *v, RGB (*color_of)(float, void *), void *ud) {
    cairo_save(cr);
    cairo_rectangle(cr, v->panel.x, v->panel.y, v->panel.w, v->panel.h);
    cairo_clip(cr);
    Box *boxes = screen_boxes(g, v);
    fill_cells(cr, g, boxes, color_of, ud);
    free(boxes);
    cairo_restore(cr);
}

typedef struct {
    double lai_2024, lai_2026, change;
    int pixels, down;
} Block;

typedef struct {
    int valid, burned, burned_no_lai;
    Block inside, outside;
} Summary;

static Block block_stats(const float *before, const float *after, const unsigned char *mask, int n, int want) {
    double sb = 0, sa = 0, sd = 0;
    int pixels = 0, down = 0;
    for (int i = 0; i < n; i++) {
        if ((mask[i] != 0) != (want != 0)) continue;
        if (!finite_f(before[i]) || !finite_f(after[i])) continue;
        double d = (double)after[i] - (double)before[i];
        sb += before[i];
        sa += after[i];
        sd += d;
        pixels++;
        if (d <= -0.3) down++;
    }
    Block b = {0};
    b.pixels = pixels;
    if (pixels) {
        b.lai_2024 = sb / pixels;
        b.lai_2026 = sa / pixels;
        b.change = sd / pixels;
        b.down = down;
    }
    return b;
}

static Summary summarize(const Grid *before, const Grid *after, const Grid *fire) {
    int n = before->w * before->h;
    float *b = before->v;
    float *a = xmalloc((size_t)n * sizeof(float));
    unsigned char *mask = xmalloc((size_t)n);
    int burned = 0, burned_no = 0, valid = 0;
    for (int r = 0; r < before->h; r++) {
        for (int c = 0; c < before->w; c++) {
            double x, y;
            cell_center(before, c, r, &x, &y);
            int i = r * before->w + c;
            a[i] = sample_nn(after, x, y);
            float day = sample_nn(fire, x, y);
            mask[i] = finite_f(day) ? 1 : 0;
            if (mask[i]) {
                burned++;
                if (!finite_f(b[i])) burned_no++;
            }
            if (finite_f(b[i]) && finite_f(a[i])) valid++;
        }
    }
    Summary s = {0};
    s.valid = valid;
    s.burned = burned;
    s.burned_no_lai = burned_no;
    s.inside = block_stats(b, a, mask, n, 1);
    s.outside = block_stats(b, a, mask, n, 0);
    free(a);
    free(mask);
    return s;
}

static Grid *aligned_float(const Grid *dst_geom, const Grid *src) {
    Grid *g = calloc(1, sizeof(Grid));
    if (!g) die("out of memory");
    *g = *dst_geom;
    g->v = xmalloc((size_t)g->w * g->h * sizeof(float));
    g->to_merc = dst_geom->to_merc;
    for (int r = 0; r < g->h; r++) {
        for (int c = 0; c < g->w; c++) {
            double x, y;
            cell_center(dst_geom, c, r, &x, &y);
            g->v[(size_t)r * g->w + c] = sample_nn(src, x, y);
        }
    }
    return g;
}

static void finish_png(cairo_surface_t *surface, const char *path) {
    if (cairo_surface_write_to_png(surface, path) != CAIRO_STATUS_SUCCESS) die("could not write PNG");
    fprintf(stderr, "wrote %s (%d x %d)\n", path, cairo_image_surface_get_width(surface),
            cairo_image_surface_get_height(surface));
}

static View view_for(Rect panel, double left, double bottom, double right, double top) {
    View v = {left, right, bottom, top, panel};
    return v;
}

static void render_burn(CURL *curl, Grid *panels[3], double left, double bottom, double right, double top,
                        cairo_font_face_t *regular, cairo_font_face_t *bold) {
    int panel_h = (int)lround(PANEL_W * (top - bottom) / (right - left));
    double gutter = 320;
    double head = 72;
    double block_title = 96;
    double vgap = 48;
    double caption_h = 80;
    int fig_w = MARGIN + PANEL_W + (int)gutter + MARGIN;
    int fig_h = MARGIN + (int)head + 3 * (int)(block_title + panel_h) + 2 * (int)vgap + (int)caption_h + MARGIN;
    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, fig_w, fig_h);
    cairo_t *cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);

    cairo_surface_t *base = fetch_basemap(curl, left, bottom, right, top, PANEL_W, panel_h);
    const char *titles[][2] = {
        {"December 2024", "no burned pixels"},
        {"January 2025, Eaton Fire", "burn day on the Altadena map"},
        {"February 2025", "no burned pixels"},
    };
    RGB black = {0.08, 0.08, 0.08};
    double map_cx = MARGIN + PANEL_W / 2.0;
    draw_centered(cr, bold, 40, map_cx, MARGIN + 42,
                  "MODIS burned area (MCD64A1, 500 m) on a map of Altadena", black);
    double y = MARGIN + head;
    Rect drawn[3];
    for (int i = 0; i < 3; i++) {
        Rect panel = {MARGIN, y + block_title, PANEL_W, panel_h};
        drawn[i] = panel;
        View v = view_for(panel, left, bottom, right, top);
        paint_basemap(cr, base, panel);
        draw_raster(cr, panels[i], &v, color_burn, NULL);
        draw_places(cr, regular, &v);
        draw_scale(cr, regular, &v);
        draw_panel_frame(cr, panel);
        draw_centered(cr, bold, 36, map_cx, panel.y - 58, titles[i][0], black);
        draw_centered(cr, regular, 28, map_cx, panel.y - 20, titles[i][1], black);
        y = panel.y + panel.h + vgap;
    }
    Rect bar = {MARGIN + PANEL_W + 48, drawn[1].y + drawn[1].h * 0.14, 36, drawn[1].h * 0.72};
    colorbar_burn(cr, regular, bar);
    draw_centered(cr, regular, 22, map_cx, fig_h - MARGIN - 16,
                  "Basemap: Esri World Topographic Map. Colored cells are 500 m burned-land pixels, not individual trees.",
                  black);
    finish_png(surface, "output/altadena_mcd64a1_burn.png");
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    cairo_surface_destroy(base);
}

static void render_lai(CURL *curl, Grid *before, Grid *after, Grid *change, const unsigned char *mask,
                       double left, double bottom, double right, double top, const Summary *summary,
                       cairo_font_face_t *regular, cairo_font_face_t *bold) {
    int panel_h = (int)lround(PANEL_W * (top - bottom) / (right - left));
    double gutter = 520;
    double head = 72;
    double block_title = 96;
    double vgap = 48;
    double caption_h = 120;
    int fig_w = MARGIN + PANEL_W + (int)gutter + MARGIN;
    int fig_h = MARGIN + (int)head + 3 * (int)(block_title + panel_h) + 2 * (int)vgap + (int)caption_h + MARGIN;
    cairo_surface_t *surface = cairo_image_surface_create(CAIRO_FORMAT_ARGB32, fig_w, fig_h);
    cairo_t *cr = cairo_create(surface);
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_paint(cr);
    cairo_surface_t *base = fetch_basemap(curl, left, bottom, right, top, PANEL_W, panel_h);

    Grid *grids[3] = {before, after, change};
    RGB (*colors[3])(float, void *) = {color_lai, color_lai, color_change};
    const char *t1[] = {"7-14 October 2024", "8-15 October 2026", "Change, 2026 minus 2024"};
    const char *t2[] = {"leaf area before the fire", "leaf area, same 8-day period", "red is canopy loss"};
    RGB black = {0.08, 0.08, 0.08};
    double map_cx = MARGIN + PANEL_W / 2.0;
    draw_centered(cr, bold, 36, map_cx, MARGIN + 42,
                  "MODIS leaf area (MYD15A2H Lai_500m, 500 m) on a map of Altadena", black);
    double y = MARGIN + head;
    Rect drawn[3];
    for (int i = 0; i < 3; i++) {
        Rect panel = {MARGIN, y + block_title, PANEL_W, panel_h};
        drawn[i] = panel;
        View v = view_for(panel, left, bottom, right, top);
        paint_basemap(cr, base, panel);
        draw_raster(cr, grids[i], &v, colors[i], NULL);
        if (i == 2) draw_outline(cr, before, &v, mask);
        draw_places(cr, regular, &v);
        draw_scale(cr, regular, &v);
        draw_panel_frame(cr, panel);
        draw_centered(cr, bold, 36, map_cx, panel.y - 58, t1[i], black);
        draw_centered(cr, regular, 28, map_cx, panel.y - 20, t2[i], black);
        y = panel.y + panel.h + vgap;
    }
    double bar_y = drawn[1].y + drawn[1].h * 0.12;
    double bar_h = drawn[1].h * 0.76;
    Rect bar1 = {MARGIN + PANEL_W + 40, bar_y, 34, bar_h};
    Rect bar2 = {bar1.x + 240, bar_y, 34, bar_h};
    colorbar_continuous(cr, regular, bar1, YLGN, 9, 0, 2.5, "Leaf area index", "%.1f");
    colorbar_continuous(cr, regular, bar2, RDYLGN, 11, -1.5, 1.5, "Change", "%.1f");

    char line1[] = "Open areas are built-up land, where this product stores 250 and does not report leaf area.";
    char line2[] = "The black line is the January 2025 burn.";
    char line3[240];
    snprintf(line3, sizeof line3,
             "Inside the burn, leaf area %.2f -> %.2f. Outside the burn, %.2f -> %.2f. Basemap: Esri World Topographic Map.",
             summary->inside.lai_2024, summary->inside.lai_2026, summary->outside.lai_2024, summary->outside.lai_2026);
    draw_centered(cr, regular, 22, map_cx, fig_h - MARGIN - 64, line1, black);
    draw_centered(cr, regular, 22, map_cx, fig_h - MARGIN - 36, line2, black);
    draw_centered(cr, regular, 22, map_cx, fig_h - MARGIN - 8, line3, black);

    finish_png(surface, "output/altadena_lai.png");
    cairo_destroy(cr);
    cairo_surface_destroy(surface);
    cairo_surface_destroy(base);
}

static void write_summary(const Summary *s) {
    FILE *f = fopen("output/altadena_lai_summary.json", "w");
    if (!f) die("could not write summary json");
    fprintf(f,
            "{\n"
            "  \"valid_leaf_area_pixels\": %d,\n"
            "  \"burned_pixels\": %d,\n"
            "  \"burned_pixels_without_lai\": %d,\n"
            "  \"inside_burn\": {\n"
            "    \"pixels\": %d,\n"
            "    \"lai_2024\": %.3f,\n"
            "    \"lai_2026\": %.3f,\n"
            "    \"change\": %.3f,\n"
            "    \"pixels_down_by_0_3\": %d\n"
            "  },\n"
            "  \"outside_burn\": {\n"
            "    \"pixels\": %d,\n"
            "    \"lai_2024\": %.3f,\n"
            "    \"lai_2026\": %.3f,\n"
            "    \"change\": %.3f,\n"
            "    \"pixels_down_by_0_3\": %d\n"
            "  }\n"
            "}\n",
            s->valid, s->burned, s->burned_no_lai, s->inside.pixels, s->inside.lai_2024, s->inside.lai_2026,
            s->inside.change, s->inside.down, s->outside.pixels, s->outside.lai_2024, s->outside.lai_2026,
            s->outside.change, s->outside.down);
    fclose(f);
}

static unsigned char *burn_mask(const Grid *before, const Grid *fire, int *burned_out) {
    int n = before->w * before->h;
    unsigned char *mask = xmalloc((size_t)n);
    int burned = 0;
    for (int r = 0; r < before->h; r++) {
        for (int c = 0; c < before->w; c++) {
            double x, y;
            cell_center(before, c, r, &x, &y);
            int on = finite_f(sample_nn(fire, x, y));
            mask[(size_t)r * before->w + c] = (unsigned char)on;
            burned += on;
        }
    }
    *burned_out = burned;
    return mask;
}

int main(void) {
    mkdir("cache", 0755);
    mkdir("output", 0755);
    curl_global_init(CURL_GLOBAL_DEFAULT);
    CURL *curl = curl_easy_init();
    if (!curl) die("curl init failed");
    g_proj = proj_context_create();
    if (!g_proj) die("proj init failed");

    FT_Library ft;
    if (FT_Init_FreeType(&ft)) die("freetype init failed");
    FT_Face regular_face = NULL, bold_face = NULL;
    const char *regular_path = "/System/Library/Fonts/Supplemental/Arial.ttf";
    const char *bold_path = "/System/Library/Fonts/Supplemental/Arial Bold.ttf";
    if (FT_New_Face(ft, regular_path, 0, &regular_face) || FT_New_Face(ft, bold_path, 0, &bold_face))
        die("could not load Arial");
    cairo_font_face_t *regular = cairo_ft_font_face_create_for_ft_face(regular_face, 0);
    cairo_font_face_t *bold = cairo_ft_font_face_create_for_ft_face(bold_face, 0);

    const char *burn_items[] = {
        "MCD64A1.A2024336.h08v05.061.2025035195913",
        "MCD64A1.A2025001.h08v05.061.2025070160807",
        "MCD64A1.A2025032.h08v05.061.2025098201953",
    };
    const char *burn_names[] = {"before", "fire", "after"};
    Grid *burn[3];
    for (int i = 0; i < 3; i++) {
        fprintf(stderr, "burn %s\n", burn_names[i]);
        burn[i] = load_grid(curl, "modis-64A1-061", burn_items[i], "Burn_Date", 0);
        int n = 0;
        for (int k = 0; k < burn[i]->w * burn[i]->h; k++)
            if (finite_f(burn[i]->v[k])) n++;
        fprintf(stderr, "  burned pixels %d\n", n);
    }

    fprintf(stderr, "lai before\n");
    Grid *before = load_grid(curl, "modis-15A2H-061", "MYD15A2H.A2024281.h08v05.061.2024290063654", "Lai_500m", 1);
    fprintf(stderr, "lai after\n");
    Grid *after = load_grid(curl, "modis-15A2H-061", "MYD15A2H.A2026281.h08v05.061.2026290035710", "Lai_500m", 1);

    Summary summary = summarize(before, after, burn[1]);
    fprintf(stderr, "inside burn %.3f -> %.3f (%d px), outside %.3f -> %.3f (%d px)\n", summary.inside.lai_2024,
            summary.inside.lai_2026, summary.inside.pixels, summary.outside.lai_2024, summary.outside.lai_2026,
            summary.outside.pixels);
    write_summary(&summary);

    double left, bottom, right, top;
    lonlat_merc(WEST, SOUTH, &left, &bottom);
    lonlat_merc(EAST, NORTH, &right, &top);

    render_burn(curl, burn, left, bottom, right, top, regular, bold);

    Grid *after_on = aligned_float(before, after);
    Grid *change = calloc(1, sizeof(Grid));
    if (!change) die("out of memory");
    *change = *before;
    change->v = xmalloc((size_t)before->w * before->h * sizeof(float));
    change->to_merc = before->to_merc;
    for (int i = 0; i < before->w * before->h; i++) {
        if (finite_f(before->v[i]) && finite_f(after_on->v[i])) change->v[i] = after_on->v[i] - before->v[i];
        else change->v[i] = NAN;
    }
    int burned = 0;
    unsigned char *mask = burn_mask(before, burn[1], &burned);
    render_lai(curl, before, after_on, change, mask, left, bottom, right, top, &summary, regular, bold);

    free(mask);
    free(change->v);
    free(change);
    free(after_on->v);
    free(after_on);
    free_grid(before);
    free_grid(after);
    for (int i = 0; i < 3; i++) free_grid(burn[i]);
    cairo_font_face_destroy(regular);
    cairo_font_face_destroy(bold);
    FT_Done_Face(regular_face);
    FT_Done_Face(bold_face);
    FT_Done_FreeType(ft);
    proj_context_destroy(g_proj);
    curl_easy_cleanup(curl);
    curl_global_cleanup();
    free(g_token);
    return 0;
}
