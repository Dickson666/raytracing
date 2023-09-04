#include<cstdio>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<limits>
#define _USE_MATH_DEFINES
#include"geometry.h"

// #define double float

using namespace std;

vec<3> bg_col;

struct Material{
    vec<3> albedo;
    vec<3> diffuse_color;
    double specular_exponent;
    Material(const vec<3>& a, const vec<3>& d, const double& s) : albedo(a), diffuse_color(d), specular_exponent(s) {}
    Material() : albedo(), diffuse_color(), specular_exponent() {}
};

vec<3> reflect(const vec<3>& L, const vec<3>& N){
    return 2 * (L * N) * N - L;
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

bool scene_intersect(const vec<3>& orig, const vec<3> &dir, vector<Sphere> spheres, vec<3>& N, vec<3>& point, int& idx){//是否有相交
    double ans = numeric_limits<double>::max();
    double ans1 = ans;
    bool res = 0;
    vec<3> N1;
    for(int i = 0; i < spheres.size(); i++)
        if(spheres[i].chk_ray(orig, dir, ans1, N1) && ans1 < ans)
            ans = ans1, N = N1, idx = i, point = spheres[i].center + N1 * spheres[i].radius, res = 1;
    return res;
}

double get_diffuse_light(const vec<3>& N, const vec<3>& point, vector<Sphere> spheres){//漫反射
    double diffuse_light_intensity = 0;
    for(int i = 0; i < lights.size(); i++){
        vec<3> light_dir = (lights[i].position - point).normalize();
        double light_distance = (lights[i].position - point).norm();
        vec<3> shadow_orig = light_dir * N < 0 ? point - (N * 1e-9) : point + (N * 1e-9);
        vec<3> shadow_N, shadow_point;
        int idx = 0;
        if(scene_intersect(shadow_orig, light_dir, spheres, shadow_N, shadow_point, idx) && (lights[i].position - shadow_point).norm() < light_distance)
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
        int idx = 0;
        if(scene_intersect(shadow_orig, light_dir, spheres, shadow_N, shadow_point, idx) && (lights[i].position - shadow_point).norm() < light_distance)
            continue;
        specular_light_intensity += lights[i].intensity * max(0.0, pow(reflect(light_dir, N) * V, alpha));
    }
    return specular_light_intensity;
}

vec<3> get_color(const vec<3>& orig, const vec<3>& dir, vector<Sphere> spheres, int dep){
    vec<3> res = bg_col;
    vec<3> N, point;
    int i = 0;
    if(dep > 4 || !scene_intersect(orig, dir, spheres, N, point, i))
        return bg_col;
    Sphere sphere = spheres[i];
    vec<3> reflect_dir = reflect(-dir, N).normalize();
    vec<3> reflect_orig = dir * N > 0 ? point - N * 1e-3 : point + N * 1e-3;
    vec<3> reflect_color = get_color(reflect_orig, reflect_dir, spheres, dep + 1);
    res = sphere.material.diffuse_color * get_diffuse_light(N, point, spheres) * sphere.material.albedo[0]
        + vec<3>{0, 255, 255, 255} * get_specular_light(N, point, (vec<3>{0} - sphere.center).normalize(), sphere.material.specular_exponent, spheres) * sphere.material.albedo[1]
        + reflect_color * sphere.material.albedo[2];
    res[1] = floor(res[1]), res[2] = floor(res[2]), res[3] = floor(res[3]);
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

    for(int i = 1; i <= width; i++)
        for(int j = 1; j <= height; j++){
            double x = (2 * (i - 0.5) / (1.0 * width) - 1) * tan(fov / 2.) * width / (double) height;
            double y = -(2 * (j - 0.5) / (1.0 * height) - 1) * tan(fov / 2.);

            vec<3> Now;
            Now[1] = x, Now[2] = y, Now[3] = -1;
            Now = Now.normalize();
            Image[i + (j - 1) * width] = (get_color(Camera, Now, spheres, 0));
        }

    ofstream ofs;
    ofs.open("./out.ppm");
    ofs << "P3\n" << width <<" "<< height<<"\n255\n";
    for(int i = 1; i <= width * height; i++){
        vec<3> c = Image[i];
        double mx = max(c[1], max(c[2], c[3]));
        if(mx > 255)
            c = c / mx * 255;
        for(int j = 1; j <= 3; j++){
            ofs << (int(c[j])) <<" ";
        }
    }
    ofs.close();
}
signed main(){
    bg_col[1] = bg_col[2] = bg_col[3] = 255;
    vec<3> ct = vec<3>{0, 0, 0, -10};
    Material m1(vec<3>{0.6, 0.3, 0.1}, vec<3>{0, 124, 100, 32}, 50);
    Material m2(vec<3>{0.9, 0.1, 0}, vec<3>{0, 5, 114, 107} , 10);
    Material m3(vec<3>{0.0, 10.0, 0.8}, vec<3>{0, 255, 255, 255}, 1425);
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
    spheres.push_back(Sphere(vec<3>{0, -1, -1.5, -12}, 2, m3));
    spheres.push_back(Sphere(vec<3>{0, 1.5, -0.5, -18}, 3, m2));
    spheres.push_back(Sphere(vec<3>{0, 7, 5, -18}, 4, m3));
    work(spheres);
}