#include<cstdio>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<limits>
#define _USE_MATH_DEFINES
#include"geometry.h"
#define STB_IMAGE_WRITE_IMPLEMENTATION
#include"stb/stb_image_write.h"
#define STB_IMAGE_IMPLEMENTATION
#include"stb/stb_image.h"

// #define double float

using namespace std;

// vec<3> bg_col;

vector<vec<3> > bgc;

int bg_width, bg_height;


struct Material{
    vec<4> albedo;
    vec<3> diffuse_color;
    double specular_exponent;
    double refractive_index;
    Material(const vec<4>& a, const vec<3>& d, const double& s, const double& r) : albedo(a), diffuse_color(d), specular_exponent(s), refractive_index(r) {}
    Material() : albedo(vec<4>{1,0,0,0}), diffuse_color(), specular_exponent(), refractive_index() {}
};

vec<3> reflect(const vec<3>& L, const vec<3>& N){
    return 2 * (L * N) * N - L;
}

vec<3> refract(vec<3> L, vec<3> N, double refractive_index){
    L = L.normalize();
    double cos1 = L * N;
    if(cos1 < 0)
        N = -N, refractive_index = 1.0 / refractive_index, cos1 = -cos1;
    double cos2 = sqrt(1 - refractive_index * refractive_index * (1 - cos1 * cos1));
    if(cos2 < 0)
        return vec<3>{0, 0, 0, 0};
    double tan1 = sqrt((1 - cos1 * cos1)) / cos1;
    double tan2 = sqrt((1 - cos2 * cos2)) / cos2;
    return (L - N * cos1) / tan1 * tan2 + N * cos1;
}

struct Sphere{
    vec<3> center;
    double radius;
    Material material;
    // vec<3> color;

    Sphere(const vec<3>& c, const double& r, const Material& m) : center(c), radius(r), material(m) {}
    Sphere() : center(), radius(), material() {}

    bool chk_ray(const vec<3>& orig, const vec<3>& dir, double& t0, vec<3>& N){// t0 是 orig 到圆的距离
        vec<3> L = center - orig;
        double dis = L * dir;
        double ctr = L * L - dis * dis;
        if(ctr > radius * radius) return false;
        double ito = sqrt(radius * radius - ctr);
        t0 = dis - ito;
        double t1 = dis + ito;
        if(t0 < 0) t0 = t1;
        if(t0 < 0) return false;
        N = ((dir * t0) - L).normalize();
        return true;
    }
};

struct Light{
    vec<3> position;
    double intensity;

    Light(const vec<3>& p, const double& i) : position(p), intensity(i){}
};

vector<Light> lights;

vec<3> get_bg_col(vec<3> dir){
    if(dir.norm() != 1)
        return vec<3>{0, 0, 0, 0};
    int x = -1, y = -1;
    y = bg_height / 2 + (dir[2] > 0 ? -1 : 1) * asin(fabs(dir[2])) / (M_PI / 2) * (bg_height / 2);
    x = (dir[1] < 0 ? 3 : 1) * bg_width / 4 + (dir[3] * dir[1] > 0 ? -1 : 1) * (bg_width / 4 * atan(fabs(dir[1]) / fabs(dir[3])) / (M_PI / 2));
    cout << dir <<" ";
    cout << x <<" " <<y <<"\n";
    return bgc[x - 1 + (y - 1) * bg_width];
}

bool scene_intersect(const vec<3>& orig, const vec<3> &dir, vector<Sphere> spheres, vec<3>& N, vec<3>& point, Material& material){//是否有相交
    double ans = numeric_limits<double>::max();
    double ans1 = ans;
    bool res = 0;
    vec<3> N1;
    for(int i = 0; i < spheres.size(); i++)
        if(spheres[i].chk_ray(orig, dir, ans1, N1) && ans1 < ans)
            ans = ans1, N = N1, material = spheres[i].material, point = spheres[i].center + N1 * spheres[i].radius, res = 1;
    if(fabs(dir[2]) <= 1e-3)
        return res;
    double dis = -(orig[2] + 4) / dir[2];
    vec<3> point1 = orig + dis * dir;
    if(dis > 0 && point1[1] >= -10 && point1[1] <= 10 && point1[3] >= -30 &&point1[3] <= -10&& dis < ans){
        res = 1;
        N = vec<3>{0, 0, 1, 0};
        point = point1;
        material.diffuse_color = (point1[1] < 0) ^ ((int(point1[1] / 2) + int(point1[3] / 2)) & 1) ? vec<3>{0, 123, 114, 9} : vec<3>{0, 9, 114, 123}; 
    }
    return res;
}

double get_diffuse_light(const vec<3>& N, const vec<3>& point, vector<Sphere> spheres){//漫反射
    double diffuse_light_intensity = 0;
    for(int i = 0; i < lights.size(); i++){
        vec<3> light_dir = (lights[i].position - point).normalize();
        double light_distance = (lights[i].position - point).norm();
        vec<3> shadow_orig = light_dir * N < 0 ? point - (N * 1e-9) : point + (N * 1e-9);
        vec<3> shadow_N, shadow_point;
        Material material;
        if(scene_intersect(shadow_orig, light_dir, spheres, shadow_N, shadow_point, material) && (lights[i].position - shadow_point).norm() < light_distance)
            continue;
        // cout << light_dir * N <<"\n";
        diffuse_light_intensity += lights[i].intensity * max(0.0, light_dir * N);
    }
    return diffuse_light_intensity;
}

double get_specular_light(const vec<3>& N, const vec<3>& point, const vec<3>& V, const double& alpha, vector<Sphere> spheres){//镜面反射
    double specular_light_intensity = 0;
    for(int i = 0; i < lights.size(); i++){
        vec<3> light_dir = (lights[i].position - point).normalize();
        double light_distance = (lights[i].position - point).norm();
        vec<3> shadow_orig = light_dir * N < 0 ? point - (N * 1e-9) : point + (N * 1e-9);
        vec<3> shadow_N, shadow_point;
        Material material;
        if(scene_intersect(shadow_orig, light_dir, spheres, shadow_N, shadow_point, material) && (lights[i].position - shadow_point).norm() < light_distance)
            continue;
        specular_light_intensity += lights[i].intensity * max(0.0, pow(reflect(light_dir, N) * V, alpha));
    }
    return specular_light_intensity;
}

vec<3> get_color(const vec<3>& orig, const vec<3>& dir, vector<Sphere> spheres, int dep){
    vec<3> res = get_bg_col(dir);
    vec<3> N, point;
    Material material;
    if(dep > 4 || !scene_intersect(orig, dir, spheres, N, point, material))
        return res;
    // Sphere sphere = spheres[i];
    vec<3> reflect_dir = reflect(-dir, N).normalize();
    vec<3> refract_dir = refract(dir, N, material.refractive_index).normalize();
    // if(sphere.material.albedo[3])
    //     cout << dir <<" " << N <<" "<<refract_dir<<"\n",system("pause");
    vec<3> reflect_orig = reflect_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    vec<3> refract_orig = refract_dir * N < 0 ? point - N * 1e-3 : point + N * 1e-3;
    vec<3> reflect_color = get_color(reflect_orig, reflect_dir, spheres, dep + 1);
    vec<3> refract_color = get_color(refract_orig, refract_dir, spheres, dep + 1);
    // if(sphere.material.albedo[3])
    //     cout << refract_color <<" dir:"<<dir<<" refract_orig:"<< refract_orig<<" refract_dir:"<< refract_dir<<" point"<< point<<"\n", system("pause");
    res = material.diffuse_color * get_diffuse_light(N, point, spheres) * material.albedo[0]
        + vec<3>{0, 255, 255, 255} * get_specular_light(N, point, (vec<3>{0} - point).normalize(), material.specular_exponent, spheres) * material.albedo[1]
        + reflect_color * material.albedo[2] + refract_color * material.albedo[3];
    res[1] = floor(res[1]), res[2] = floor(res[2]), res[3] = floor(res[3]);
    // if(sphere.material.albedo[3] && res != bg_col && res != vec<3>{0, 204, 204, 204})
    //     cout << res <<"\n"/*, system("pause")*/;
    return res;
}

// Sphere sphere1;

void work(vector<Sphere> spheres){
    int height = 1080;
    int width = 1920;
    double fov = M_PI / 3.0;

    vec<3> Camera;

    Camera[1] = 0, Camera[2] = 0, Camera[3] = 0;

    vector<vec<3> >Image(width * height + 1);

    for(int i = 1; i <= width; i++){
        cout << i <<"\n";
        for(int j = 1; j <= height; j++){
            double x = (2 * (i - 0.5) / (1.0 * width) - 1) * tan(fov / 2.) * width / (double) height;
            double y = -(2 * (j - 0.5) / (1.0 * height) - 1) * tan(fov / 2.);

            vec<3> Now;
            Now[1] = x, Now[2] = y, Now[3] = -1;
            Now = Now.normalize();
            Image[i + (j - 1) * width] = (get_color(Camera, Now, spheres, 0));
        }
    }

    // ofstream ofs;
    // ofs.open("./out.ppm");
    // ofs << "P3\n" << width <<" "<< height<<"\n255\n";
    vector<unsigned char> pixmap(width * height *3);
    for(int i = 1; i <= width * height; i++){
        vec<3> c = Image[i];
        double mx = max(c[1], max(c[2], c[3]));
        if(mx > 255)
            c = c / mx * 255;
        for(int j = 1; j <= 3; j++)
            pixmap[(i - 1) * 3 + j - 1] = (unsigned char)c[j];
        // for(int j = 1; j <= 3; j++){
        //     ofs << (int(c[j])) <<" ";
        // }
    }
    stbi_write_jpg("out.jpg", width, height, 3, pixmap.data(), 100);
    // ofs.close();
}
signed main(){
    int Type = -1;
    unsigned char *background = stbi_load("background1.jpg", &bg_width, &bg_height, &Type, 0);
    if(!background || Type != 3){
        printf("No image");
        return -1;
    }
    vector<vec<3> >bg_col(bg_width * bg_height);
    for(int i = 0; i < bg_width * bg_height; i += 3)
        bg_col[i] = vec<3>{0, 1.0 * background[i], 1.0 * background[i + 1], 1.0 * background[i + 2]};

    bgc = bg_col;
    // bg_col[1] = bg_col[2] = bg_col[3] = 255;
    vec<3> ct = vec<3>{0, 0, 0, -10};
    Material m1(vec<4>{0.6, 0.3, 0.1, 0}, vec<3>{0, 124, 100, 32}, 50, 1);
    Material m2(vec<4>{0.9, 0.1, 0, 0}, vec<3>{0, 5, 114, 107} , 10, 1);
    Material m3(vec<4>{0.0, 10.0, 0.8, 0}, vec<3>{0, 255, 255, 255}, 1425, 1);
    Material m4(vec<4>{0.0, 0.5, 0.1, 0.8}, vec<3>{0, 135, 175, 200}, 125, 1.5);
    // vec<3> co1 = vec<3>{0, 124, 100, 0};
    lights.push_back(Light(vec<3>{0, -20, 20, 20}, 0.7));
    lights.push_back(Light(vec<3>{0, 30, 50, -25}, 1.3));
    lights.push_back(Light(vec<3>{0, 30, 20, 30},  1.2));
    // co1[1] = 76, co1[2] = 208, co1[3] = 212;
    // ct[1] = 0, ct[2] = 0, ct[3] = -10;
    vector<Sphere> spheres;
    // vector<Material> materials;
    // spheres.push_back(Sphere(ct, 6, co1));
    spheres.push_back(Sphere(vec<3>{0, -3, 0, -16}, 2, m1));
    spheres.push_back(Sphere(vec<3>{0, -1, -1.5, -12}, 2, m4));
    spheres.push_back(Sphere(vec<3>{0, 1.5, -0.5, -18}, 3, m2));
    spheres.push_back(Sphere(vec<3>{0, 7, 5, -18}, 4, m3));
    work(spheres);
}