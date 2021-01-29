/*
 * A simple test program to display the output of the voronoi generator

VERSION
    0.2     2017-04-16  - Added support for reading .csv files
    0.1                 - Initial version

 */

#include <stdlib.h>
#include <stdint.h>
#include <stdio.h> // printf
#include <stdbool.h>
#include <iostream>
#include <string>
#include <functional>
#include <memory>
#include <vector>
#include <unordered_map>
#include <algorithm>

#if defined(_MSC_VER)
#include <malloc.h>
#define alloca _alloca
#else
#include <alloca.h>
#endif

#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"

int wrap_stbi_write_png(char const *filename, int w, int h, int comp, const void *data, int stride_in_bytes)
{
    return stbi_write_png(filename, w, h, comp, data, stride_in_bytes);
}


#define JC_VORONOI_IMPLEMENTATION
// If you wish to use doubles
//#define JCV_REAL_TYPE double
//#define JCV_FABS fabs
//#define JCV_ATAN2 atan2
#include "jc_voronoi.h"

#define JC_VORONOI_CLIP_IMPLEMENTATION
#include "jc_voronoi_clip.h"

// 加入hash map用于美化边界效果
#include "hashmap.h"

#ifdef HAS_MODE_FASTJET
#include <vector>
#include "../test/fastjet/voronoi.h"
#endif


static int num_pts = 10;

static void plot(int x, int y, unsigned char* image, int width, int height, int nchannels, unsigned char* color)
{
    if( x < 0 || y < 0 || x > (width-1) || y > (height-1) )
        return;
    int index = y * width * nchannels + x * nchannels;
    for( int i = 0; i < nchannels; ++i )
    {
        image[index+i] = color[i];
    }
}

// http://members.chello.at/~easyfilter/bresenham.html
static void draw_line(int x0, int y0, int x1, int y1, unsigned char* image, int width, int height, int nchannels, unsigned char* color)
{
    int dx =  abs(x1-x0), sx = x0<x1 ? 1 : -1;
    int dy = -abs(y1-y0), sy = y0<y1 ? 1 : -1;
    int err = dx+dy, e2; // error value e_xy

    for(;;)
    {  // loop
        plot(x0,y0, image, width, height, nchannels, color);
        if (x0==x1 && y0==y1) break;
        e2 = 2*err;
        if (e2 >= dy) { err += dy; x0 += sx; } // e_xy+e_x > 0
        if (e2 <= dx) { err += dx; y0 += sy; } // e_xy+e_y < 0
    }
}

// http://fgiesen.wordpress.com/2013/02/08/triangle-rasterization-in-practice/
static inline int orient2d(const jcv_point* a, const jcv_point* b, const jcv_point* c)
{
    return ((int)b->x - (int)a->x)*((int)c->y - (int)a->y) - ((int)b->y - (int)a->y)*((int)c->x - (int)a->x);
}

static inline int min2(int a, int b)
{
    return (a < b) ? a : b;
}

static inline int max2(int a, int b)
{
    return (a > b) ? a : b;
}

static inline int min3(int a, int b, int c)
{
    return min2(a, min2(b, c));
}
static inline int max3(int a, int b, int c)
{
    return max2(a, max2(b, c));
}

static void draw_triangle(const jcv_point* v0, const jcv_point* v1, const jcv_point* v2, unsigned char* image, int width, int height, int nchannels, unsigned char* color)
{
    int area = orient2d(v0, v1, v2);
    if( area == 0 )
        return;

    // Compute triangle bounding box
    int minX = min3((int)v0->x, (int)v1->x, (int)v2->x);
    int minY = min3((int)v0->y, (int)v1->y, (int)v2->y);
    int maxX = max3((int)v0->x, (int)v1->x, (int)v2->x);
    int maxY = max3((int)v0->y, (int)v1->y, (int)v2->y);

    // Clip against screen bounds
    minX = max2(minX, 0);
    minY = max2(minY, 0);
    maxX = min2(maxX, width - 1);
    maxY = min2(maxY, height - 1);

    // Rasterize
    jcv_point p;
    for (p.y = (jcv_real)minY; p.y <= maxY; p.y++) {
        for (p.x = (jcv_real)minX; p.x <= maxX; p.x++) {
            // Determine barycentric coordinates
            int w0 = orient2d(v1, v2, &p);
            int w1 = orient2d(v2, v0, &p);
            int w2 = orient2d(v0, v1, &p);

            // If p is on or inside all edges, render pixel.
            if (w0 >= 0 && w1 >= 0 && w2 >= 0)
            {
                plot((int)p.x, (int)p.y, image, width, height, nchannels, color);
            }
        }
    }
}

static void relax_points(const jcv_diagram* diagram, jcv_point* points)
{
    const jcv_site* sites = jcv_diagram_get_sites(diagram);
    for( int i = 0; i < diagram->numsites; ++i )
    {
        const jcv_site* site = &sites[i];
        jcv_point sum = site->p;
        int count = 1;

        const jcv_graphedge* edge = site->edges;

        while( edge )
        {
            sum.x += edge->pos[0].x;
            sum.y += edge->pos[0].y;
            ++count;
            edge = edge->next;
        }

        points[site->index].x = sum.x / count;
        points[site->index].y = sum.y / count;
    }
}

static void Usage()
{
    printf("Usage: main [options]\n");
    printf("\t-n <num points>\n");
    printf("\t-r <num relaxations>\n");
    printf("\t-i <inputfile>\t\tA list of 2-tuples (float, float) representing 2-d coordinates\n");
    printf("\t-o <outputfile.png>\n");
    printf("\t-w <width>\n");
    printf("\t-h <height>\n");
}

// Search for any of the common characters: \n,;
static inline int is_csv(const char* chars, uint32_t len)
{
    for( uint32_t i = 0; i < len; ++i )
    {
        char c = chars[i];
        if( c == '\n' || c == ',' || c == ' ' || c == ';' || c == '\t' )
            return 1;
    }
    return 0;
}

static int debug_skip_point(const jcv_point* pt)
{
    (void)pt;
    return 0;
}

static int read_input(const char* path, jcv_point** points, uint32_t* length, jcv_rect** rect)
{
    if( !path )
    {
        return 1;
    }

    FILE* file = 0;
    if( strcmp(path, "-") == 0 )
        file = stdin;
    else
        file = fopen(path, "rb");

    if( !file )
    {
        fprintf(stderr, "Failed to open %s for reading\n", path);
        *length = 0;
        return 1;
    }

    uint32_t capacity = 0;
    uint32_t len = 0;
    jcv_point* pts = 0;

    int mode = -1;
    const uint32_t buffersize = 32627;
    char* buffer = (char*)alloca(buffersize);
    uint32_t bufferoffset = 0;

    while( !feof(file) )
    {
        size_t num_read = fread((void*)&buffer[bufferoffset], 1, buffersize - bufferoffset, file);
        num_read += bufferoffset;

        if( mode == -1 )
        {
            mode = is_csv(buffer, (uint32_t)num_read);
        }

        if( mode == 0 ) // binary
        {
            uint32_t num_points = (uint32_t) num_read / sizeof(jcv_point);
            if( capacity < (len + num_points))
            {
                capacity += 1024;
                pts = (jcv_point*)realloc(pts, sizeof(jcv_point) * capacity);
            }
            for( uint32_t i = 0; i < num_points; ++i )
            {
                jcv_point* pt = &((jcv_point*)buffer)[i];
                if( debug_skip_point(pt) )
                {
                    continue;
                }
                pts[len].x = pt->x;
                pts[len].y = pt->y;
                ++len;
            }
            bufferoffset = (uint32_t) num_read - num_points * sizeof(jcv_point);
            memmove(buffer, &buffer[num_points * sizeof(jcv_point)], bufferoffset);
        }
        else if( mode == 1 ) // CSV mode
        {
            char* p = buffer;
            char* end = &buffer[num_read];

            while( p < end )
            {
                char* r = p;
                int end_of_line = 0;
                while( r < end )
                {
                    if (*r == '\0')
                    {
                        end_of_line = 1;
                        break;
                    }
                    if (*r == '\n')
                    {
                        end_of_line = 1;
                        *r = 0;
                        r += 1;
                        break;
                    }
                    else if( (*r == '\r' && *(r+1) == '\n') )
                    {
                        end_of_line = 1;
                        *r = 0;
                        r += 2;
                        break;
                    }
                    else if( *r == ',' || *r == ';' || *r == ':' )
                    {
                        *r = ' ';
                    }
                    ++r;
                }

                if( end_of_line )
                {
                    jcv_point pt1;
                    jcv_point pt2;
                    int numscanned = sscanf(p, "%f %f %f %f\n", &pt1.x, &pt1.y, &pt2.x, &pt2.y);
                    if( numscanned == 4 )
                    {
                        if (rect)
                        {
                            *rect = static_cast<jcv_rect*>(
                                malloc(sizeof(jcv_rect)));
                            (*rect)->min = pt1;
                            (*rect)->max = pt2;
                        }
                        p = r;
                    }
                    else if( numscanned == 2 )
                    {
                        if( debug_skip_point(&pt1) )
                        {
                            continue;
                        }
                        if( capacity < (len + 1))
                        {
                            capacity += 1024;
                            pts = (jcv_point*)realloc(pts, sizeof(jcv_point) * capacity);
                        }

                        pts[len].x = pt1.x;
                        pts[len].y = pt1.y;
                        ++len;
                        p = r;
                    }
                    else
                    {
                        fprintf(stderr, "Failed to read point on line %u\n", len);
                        return 1;
                    }
                }
                else
                {
                    bufferoffset = (uint32_t)(uintptr_t) (r - p);
                    memmove(buffer, p, bufferoffset);
                    break;
                }
            }
        }
    }

    printf("Read %d points from %s\n", len, path);

    if( strcmp(path, "-") != 0 )
        fclose(file);

    *points = pts;
    *length = len;

    return 0;
}

// Remaps the point from the input space to image space
static inline jcv_point remap(const jcv_point* pt, const jcv_point* min, const jcv_point* max, const jcv_point* scale)
{
    jcv_point p;
    p.x = (pt->x - min->x)/(max->x - min->x) * scale->x;
    p.y = (pt->y - min->y)/(max->y - min->y) * scale->y;
    return p;
}

// 处理diagram 导出并优化边缘
static std::string hash_sites(const jcv_site *a, const jcv_site *b)
{
    size_t t_a = reinterpret_cast<size_t>(a);
    size_t t_b = reinterpret_cast<size_t>(b);

    if (t_a > t_b)
    {
        size_t temp_a = t_b;
        t_b = t_a;
        t_a = temp_a;
    }

    char buf[100];
    snprintf(buf, 100, "%zu%zu", t_a, t_b);
    std::string hash_str = buf;

    return hash_str;
}

static float uniform_random()
{
    float r = ((float) rand() / (RAND_MAX));
    return r;
}

// static jcv_point random_point(const jcv_site *site, const jcv_graphedge *edge, bool is_neighbor)
// {
//     jcv_point site_point;
//     if(!is_neighbor)
//     {
//         site_point = site->p;
//     } else {
//         site_point = edge->neighbor->p;
//     }

//     float epsilon = 0.02;
//     float a1 = uniform_random();
//     float a2 = uniform_random();
//     a1 = a1 > 0 ? a1 : epsilon;
//     a2 = a2 > 0 ? a2 : epsilon;

//     if (a1+a2 > 1.0)
//     {
//         a1 = 1.0 - a2;
//         a2 = 1.0 - a1;
//     }
//     // printf("generate random points\n");
//     // printf("site: %f %f\n",site_point.x, site_point.y);
//     // printf("st: %f %f\n", edge->pos[0].x, edge->pos[0].y);
//     // printf("ed: %f %f\n", edge->pos[1].x, edge->pos[1].y);
//     //printf("a1, a2: %f %f\n");
//     jcv_point v1;
//     v1.x = site_point.x - edge->pos[0].x;
//     v1.y = site_point.y - edge->pos[0].y;

//     jcv_point v2;
//     v2.x = edge->pos[1].x - edge->pos[0].x;
//     v2.y = edge->pos[1].y - edge->pos[0].y;

//     jcv_point off;
//     off.x = v1.x*a1 + v2.x*a2;
//     off.y = v1.y*a1 + v2.y*a2;

//     jcv_point result;
//     float alter_ratio = 0.4;
//     result.x = edge->pos[0].x + off.x*alter_ratio;
//     result.y = edge->pos[0].y + off.y*alter_ratio;
//     // printf("result: %f %f\n", result.x, result.y);

//     return result;
// }

static jcv_point random_point(const jcv_point st_p, jcv_point ed_p, jcv_point site_p)
{
    float epsilon = 0.2;
    float a1 = uniform_random();
    float a2 = uniform_random();
    a1 = a1 > epsilon ? a1 : epsilon;
    a2 = a2 > epsilon ? a2 : epsilon;

    if (a1+a2 > 1.0)
    {
        a1 = 1.0 - a2;
        a2 = 1.0 - a1;
    }
    // printf("generate random points\n");
    // printf("site: %f %f\n",site_point.x, site_point.y);
    // printf("st: %f %f\n", edge->pos[0].x, edge->pos[0].y);
    // printf("ed: %f %f\n", edge->pos[1].x, edge->pos[1].y);
    //printf("a1, a2: %f %f\n");
    jcv_point v1;
    v1.x = site_p.x - st_p.x;
    v1.y = site_p.y - st_p.y;

    jcv_point v2;
    v2.x = ed_p.x - st_p.x;
    v2.y = ed_p.y - st_p.y;

    jcv_point off;
    off.x = v1.x*a1 + v2.x*a2;
    off.y = v1.y*a1 + v2.y*a2;

    jcv_point result;
    float alter_ratio = 0.35;
    result.x = st_p.x + off.x*alter_ratio;
    result.y = st_p.y + off.y*alter_ratio;
    // printf("result: %f %f\n", result.x, result.y);

    return result;
}

float area2(jcv_point p, jcv_point q, jcv_point s)
{
    return p.x*q.y - p.y*q.x + q.x*s.y - q.y*s.x + s.x*p.y - s.y*p.x;
}

using alter_points = std::vector<jcv_point>;
void iter_generate_points(jcv_point p,jcv_point end_p, jcv_point site0, jcv_point site1, size_t total, alter_points &points, bool is_neighbor)
{
    // 递归升成点的函数
    if (total == 0)
        return;
    
    // 如果非凸包则不生成点
    float a0 = area2(site0, site1, p);
    float a1 = area2(site0, site1, end_p);
    if (a0*a1 >= 0.0)
        return;
    
    total = total - 1;
    
    jcv_point site_point;
    if (is_neighbor)
        site_point = site1;
    else
        site_point = site0;
    // printf("site0: %f %f\n",site0.x, site0.y);
    // printf("site1: %f %f\n",site1.x, site1.y);
    // printf("check site_point: %f %f\n", site_point.x, site_point.y);
    jcv_point new_pt = random_point(p, end_p, site_point);
    points.push_back(new_pt);
    iter_generate_points(new_pt, end_p, site0, site1, total, points, is_neighbor);
}

// 根据graphedge构造一个生成边界点的函数
static alter_points generate_alter_points(std::unordered_map<std::string, alter_points> &hash_map, const jcv_site *site, const jcv_graphedge *edge, bool &is_new)
{
    jcv_site *other_site = edge->neighbor;
    std::string hash_str = hash_sites(site, other_site);

    auto it_find = hash_map.find(hash_str);
    if(it_find != hash_map.end())
    {
        is_new = false;
        return it_find->second;
    } else {
        is_new = true;
    }
    alter_points points;
    if(edge->neighbor == nullptr)
        return points;

    jcv_point site0 = site->p;
    jcv_point site1 = edge->neighbor->p;
    size_t total = rand() % num_pts + 5;
    bool is_neighbor = (rand() % 2 > 0);
    iter_generate_points(
        edge->pos[0],
        edge->pos[1],
        site0,
        site1,
        total,
        points,
        is_neighbor
    );
    // size_t max_t = 3;
    // size_t p_cnt = rand() % max_t + 1;
    // for (size_t i=0; i<p_cnt; i++)
    // {
    //     jcv_point new_pt;
    //     bool is_neighbor = (rand() % 2 > 0);
    //     new_pt = random_point(site, edge, is_neighbor);
    //     points.push_back(new_pt);
    // }
    auto compare_func = [=](jcv_point p, jcv_point q) {
        jcv_point s = site->p;
        float area2 = p.x*q.y - p.y*q.x + q.x*s.y - q.y*s.x + s.x*p.y - s.y*p.x;
        return area2 > 0;
    };
    std::sort(points.begin(), points.end(), compare_func);
    // 将生成的points插入 hash_map
    hash_map.insert({hash_str, points});
    return points;
}
// 创建一个函数将site和edge数据存储起来
static void store_voronoi_diagram(jcv_diagram diagram, size_t width, size_t height)
{
    size_t alter_imagesize = (size_t)(width*height*3);
    unsigned char* alter_image = (unsigned char*)malloc(alter_imagesize);
    memset(alter_image, 0, alter_imagesize);

    jcv_point dimensions;
    dimensions.x = (jcv_real)width;
    dimensions.y = (jcv_real)height;

    FILE* voronoi_file = fopen("voronoi_file", "w+");
    // 首先存储site
    const jcv_site* sites = jcv_diagram_get_sites( &diagram );
    // 创建key：相邻site的hash，和vector:为新生成的点
    std::unordered_map<std::string, alter_points> hash_table;
    size_t test_cnt = 0;
    size_t site_cnt = 0;
    for ( int i = 0; i < diagram.numsites; ++i)
    {
        const jcv_site* site = &sites[i];
        // 这里先hard code 先去除包围点
        if(site->p.x >= 499.0 || site->p.x <= -499.0)
        {
            continue;
        }
        // 控制点
        fprintf(voronoi_file, "%f %f\n", site->p.x, site->p.y);
        // printf("site: %f %f\n", site->p.x, site->p.y);
        // 输入文件对应点的index
        fprintf(voronoi_file, "%d\n", site->index);
        // 边
        const jcv_graphedge* e = site->edges;
        bool is_first_edge = true;

        // 测试代码
        site_cnt += 1;
        // if(site_cnt != 5)
        //     continue;
        
        unsigned char color_tri[3];
        unsigned char basecolor = 120;
        color_tri[0] = basecolor + (unsigned char)(rand() % (235 - basecolor));
        color_tri[1] = basecolor + (unsigned char)(rand() % (235 - basecolor));
        color_tri[2] = basecolor + (unsigned char)(rand() % (235 - basecolor));

        std::vector<jcv_point> edge_pts;
        jcv_point s = remap(
            &site->p,
            &diagram.min,
            &diagram.max,
            &dimensions
        );
        while(e)
        {
            // 绘制triangle
            size_t tri_idx = 0;
            for(auto it = edge_pts.begin(); it != (edge_pts.end()-1); it++)
            {
                if (edge_pts.size() <= 0)
                    break;
                
                tri_idx += 1;
                // if (tri_idx != 2)
                //     continue;
                
                jcv_point pos0 = *it;
                jcv_point pos1 = *(it+1);
                jcv_point p0 = remap(
                    &pos0,
                    &diagram.min,
                    &diagram.max,
                    &dimensions
                );
                jcv_point p1 = remap(
                    &pos1,
                    &diagram.min,
                    &diagram.max,
                    &dimensions
                );
                // printf("in draw triangle :\n");
                // printf("site: %f %f\n", site->p.x, site->p.y);
                // printf("p0: %f %f\n", pos0.x, pos0.y);
                // printf("p1: %f %f\n", pos1.x, pos1.y);
                draw_triangle(
                    &s, &p0, &p1, alter_image, width, height, 3, color_tri
                );
            }

            edge_pts.clear();
            bool is_new = false;
            if (is_first_edge)
            {
                fprintf(voronoi_file, "%f %f\n", e->pos[0].x, e->pos[0].y);
                edge_pts.push_back(e->pos[0]);
                is_first_edge = false;
            }
            bool is_last_edge = false;
            if (e->next)
                is_last_edge = false;
            else
                is_last_edge = true;

            // printf("edge: \n");
            // printf("start: %f %f\n", e->pos[0].x, e->pos[0].y);
            // printf("end: %f %f\n", e->pos[1].x, e->pos[1].y);

            edge_pts.push_back(e->pos[0]);
            alter_points points = generate_alter_points(
                hash_table,
                site,
                e,
                is_new
            );

            if(is_new)
            {
                // 正序
                for (auto it=points.begin(); it != points.end(); it++)
                {
                    fprintf(voronoi_file, "%f %f\n", it->x, it->y);
                    edge_pts.push_back(*it);
                    // printf("insert pt: %f %f\n", it->x, it->y);
                    test_cnt++;
                }
            }
            else
            {
                // 倒序
                for (auto it=points.crbegin(); it != points.crend(); it++)
                {
                    fprintf(voronoi_file, "%f %f\n", it->x, it->y);
                    edge_pts.push_back(*it);
                    // printf("insert pt: %f %f\n", it->x, it->y);
                    test_cnt++;
                }
            }
            edge_pts.push_back(e->pos[1]);

            if (!is_last_edge)
            {   
                fprintf(voronoi_file, "%f %f\n", e->pos[1].x, e->pos[1].y);
            }
            e = e->next;
        }

        for(auto it = edge_pts.begin(); it != (edge_pts.end()-1); it++)
        {
            jcv_point pos0 = *it;
            jcv_point pos1 = *(it+1);
            jcv_point p0 = remap(
                &pos0,
                &diagram.min,
                &diagram.max,
                &dimensions
            );
            jcv_point p1 = remap(
                &pos1,
                &diagram.min,
                &diagram.max,
                &dimensions
            );
            // printf("out draw triangle :\n");
            // printf("p0: %f %f\n", pos0.x, pos0.y);
            // printf("p1: %f %f\n", pos1.x, pos1.y);
            draw_triangle(
                &s, &p0, &p1, alter_image, width, height, 3, color_tri
            );
        }
        // break;
        fprintf(voronoi_file, "\n");
    }
    printf("insert pts cnt :%zu\n", test_cnt);
    fclose(voronoi_file);

    // flip image
    int stride = width*3;
    uint8_t* row = (uint8_t*)malloc((size_t)stride);
    for( int y = 0; y < height/2; ++y )
    {
        memcpy(row, &alter_image[y*stride], (size_t)stride);
        memcpy(&alter_image[y*stride], &alter_image[(height-1-y)*stride], (size_t)stride);
        memcpy(&alter_image[(height-1-y)*stride], row, (size_t)stride);
    }

    char path[512];
    sprintf(path, "%s", "alter_img.png");
    wrap_stbi_write_png(path, width, height, 3, alter_image, stride);
    printf("wrote %s\n", path);

    free(alter_image);
}

int main(int argc, const char** argv)
{
    // Number of sites to generate
    int count = 200;
    // Image dimension
    int width = 512;
    int height = 512;
    int numrelaxations = 0;
    int mode = 0;
    const char* inputfile = 0;
    const char* clipfile = 0; // a file with clipping points
    const char* outputfile = "example.png";

    if( argc == 1 )
    {
        Usage();
        return 1;
    }

    for( int i = 1; i < argc; ++i )
    {
        if(strcmp(argv[i], "-i") == 0)
        {
            if( i+1 < argc )
                inputfile = argv[i+1];
            else
            {
                Usage();
                return 1;
            }
        }
        else if(strcmp(argv[i], "-o") == 0)
        {
            if( i+1 < argc )
                outputfile = argv[i+1];
            else
            {
                Usage();
                return 1;
            }
        }
        else if(strcmp(argv[i], "-n") == 0)
        {
            if( i+1 < argc )
                count = (int)atol(argv[i+1]);
            else
            {
                Usage();
                return 1;
            }
        }
        else if(strcmp(argv[i], "-w") == 0)
        {
            if( i+1 < argc )
                width = (int)atol(argv[i+1]);
            else
            {
                Usage();
                return 1;
            }
        }
        else if(strcmp(argv[i], "-h") == 0)
        {
            if( i+1 < argc )
                height = (int)atol(argv[i+1]);
            else
            {
                Usage();
                return 1;
            }
        }
        else if(strcmp(argv[i], "-r") == 0)
        {
            if( i+1 < argc )
                numrelaxations = (int)atol(argv[i+1]);
            else
            {
                Usage();
                return 1;
            }
        }
        else if(strcmp(argv[i], "-m") == 0)
        {
            if( i+1 < argc )
                mode = (int)atol(argv[i+1]);
            else
            {
                Usage();
                return 1;
            }
        }
        if(strcmp(argv[i], "-c") == 0)
        {
            if( i+1 < argc )
                clipfile = argv[i+1];
            else
            {
                Usage();
                return 1;
            }
        }
        else if(strcmp(argv[i], "-?") == 0 || strcmp(argv[i], "--help") == 0)
        {
            Usage();
            return 0;
        }
    }

    jcv_point* points = 0;
    jcv_rect* rect = 0;

    if( inputfile )
    {
        if( read_input(inputfile, &points, (uint32_t*)&count, &rect) )
        {
            fprintf(stderr, "Failed to read from %s\n", inputfile);
            return 1;
        }
    }
    else
    {
        points = (jcv_point*)malloc( sizeof(jcv_point) * (size_t)count);
        if( !points )
            return 1;

        int pointoffset = 10; // move the points inwards, for aestetic reasons

        srand(0);

        for( int i = 0; i < count; ++i )
        {
            points[i].x = (float)(pointoffset + rand() % (width-2*pointoffset));
            points[i].y = (float)(pointoffset + rand() % (height-2*pointoffset));
        }
    }


    jcv_point* clippoints = 0;
    int clipcount = 0;
    if( clipfile )
    {
        if( read_input(clipfile, &clippoints, (uint32_t*)&clipcount, 0) )
        {
            fprintf(stderr, "Failed to read from %s\n", clipfile);
            return 1;
        }
    }

    printf("Width/Height is %d, %d\n", width, height);
    printf("Count is %d, num relaxations is %d\n", count, numrelaxations);


    jcv_clipping_polygon polygon;
    jcv_clipper* clipper = 0;
    if (clippoints)
    {

        printf("Clip polygon '%s' used\n", clipfile);
        polygon.num_points = clipcount;
        polygon.points = clippoints;

        jcv_clipper polygonclipper;
        polygonclipper.test_fn = jcv_clip_polygon_test_point;
        polygonclipper.clip_fn = jcv_clip_polygon_clip_edge;
        polygonclipper.fill_fn = jcv_clip_polygon_fill_gaps;
        polygonclipper.ctx = &polygon;

        clipper = &polygonclipper;
    } else {
        polygon.num_points = 0;
        polygon.points = 0;
    }

    for( int i = 0; i < numrelaxations; ++i )
    {
        jcv_diagram diagram;
        memset(&diagram, 0, sizeof(jcv_diagram));
        jcv_diagram_generate(count, (const jcv_point*)points, rect, clipper, &diagram);

        relax_points(&diagram, points);

        jcv_diagram_free( &diagram );
    }

    size_t imagesize = (size_t)(width*height*3);
    unsigned char* image = (unsigned char*)malloc(imagesize);
    memset(image, 0, imagesize);

    unsigned char color_pt[] = {255, 255, 255};
    unsigned char color_line[] = {220, 220, 220};

    jcv_diagram diagram;
    jcv_point dimensions;
    dimensions.x = (jcv_real)width;
    dimensions.y = (jcv_real)height;
    {
        memset(&diagram, 0, sizeof(jcv_diagram));
        jcv_diagram_generate(count, (const jcv_point*)points, rect, clipper, &diagram);

        // If you want to draw triangles, or relax the diagram,
        // you can iterate over the sites and get all edges easily
        const jcv_site* sites = jcv_diagram_get_sites( &diagram );
        for( int i = 0; i < diagram.numsites; ++i )
        {
            const jcv_site* site = &sites[i];

            srand((unsigned int)site->index); // for generating colors for the triangles

            unsigned char color_tri[3];
            unsigned char basecolor = 120;
            color_tri[0] = basecolor + (unsigned char)(rand() % (235 - basecolor));
            color_tri[1] = basecolor + (unsigned char)(rand() % (235 - basecolor));
            color_tri[2] = basecolor + (unsigned char)(rand() % (235 - basecolor));

            jcv_point s = remap(&site->p, &diagram.min, &diagram.max, &dimensions );

            const jcv_graphedge* e = site->edges;
            while( e )
            {
                jcv_point p0 = remap(&e->pos[0], &diagram.min, &diagram.max, &dimensions );
                jcv_point p1 = remap(&e->pos[1], &diagram.min, &diagram.max, &dimensions );

                draw_triangle( &s, &p0, &p1, image, width, height, 3, color_tri);
                e = e->next;
            }
        }

        // If all you need are the edges
        const jcv_edge* edge = jcv_diagram_get_edges( &diagram );
        while( edge )
        {
            jcv_point p0 = remap(&edge->pos[0], &diagram.min, &diagram.max, &dimensions );
            jcv_point p1 = remap(&edge->pos[1], &diagram.min, &diagram.max, &dimensions );
            draw_line((int)p0.x, (int)p0.y, (int)p1.x, (int)p1.y, image, width, height, 3, color_line);
            edge = jcv_diagram_get_next_edge(edge);
        }
        store_voronoi_diagram(diagram, width, height);
        jcv_diagram_free( &diagram );
    }

    // draw the clipping polygon
    for (int i = 0; i < polygon.num_points; ++i)
    {
        jcv_point p0 = remap(&polygon.points[i], &diagram.min, &diagram.max, &dimensions );
        jcv_point p1 = remap(&polygon.points[(i+1)%polygon.num_points], &diagram.min, &diagram.max, &dimensions );
        draw_line((int)p0.x, (int)p0.y, (int)p1.x, (int)p1.y, image, width, height, 3, color_line);
    }

    // Plot the sites
    for( int i = 0; i < count; ++i )
    {
        if (clipper && !clipper->test_fn(clipper, points[i]))
            continue;
        jcv_point p = remap(&points[i], &diagram.min, &diagram.max, &dimensions );
        plot((int)p.x, (int)p.y, image, width, height, 3, color_pt);
    }

    free(clippoints);
    free(points);
    free(rect);

    // flip image
    int stride = width*3;
    uint8_t* row = (uint8_t*)malloc((size_t)stride);
    for( int y = 0; y < height/2; ++y )
    {
        memcpy(row, &image[y*stride], (size_t)stride);
        memcpy(&image[y*stride], &image[(height-1-y)*stride], (size_t)stride);
        memcpy(&image[(height-1-y)*stride], row, (size_t)stride);
    }

    char path[512];
    sprintf(path, "%s", outputfile);
    wrap_stbi_write_png(path, width, height, 3, image, stride);
    printf("wrote %s\n", path);

    free(image);

    return 0;
}
