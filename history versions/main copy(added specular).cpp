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
    vec<2> albedo;
    vec<3> diffuse_color;
    double specular_exponent;
    Material(const vec<2>& a, const vec<3>& d, const double& s) : albedo(a), diffuse_color(d), specular_exponent(s) {}
    Material() : albedo(), diffuse_color(), specular_exponent() {}
};

vec<3> Rm(const vec<3>& L, const vec<3>& N){
    return 2 * (L * N) * N - L;
}

struct Sphere{
    vec<3> center;
    double radius;
    Material material;
    // vec<3> color;

    Sphere(const vec<3>& c, const double& r, const Material& m) : center(c), radius(r), material(m) {}

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

double get_diffuse_light(const vec<3>& N, const vec<3>& point){
    double diffuse_light_intensity = 0;
    for(int i = 0; i < lights.size(); i++){
        vec<3> light_dir = (lights[i].position - point).normalize();
        // cout << light_dir * N <<"\n";
        diffuse_light_intensity += lights[i].intensity * max(0.0, light_dir * N);
    }
    return diffuse_light_intensity;
}

double get_specular_light(const vec<3>& N, const vec<3>& point, const vec<3>& V, const double& alpha){
    double specular_light_intensity = 0;
    for(int i = 0; i < lights.size(); i++){
        vec<3> light_dir = (lights[i].position - point).normalize();
        specular_light_intensity += lights[i].intensity * max(0.0, pow(Rm(light_dir, N) * V, alpha));
    }
    return specular_light_intensity;
}

vec<3> get_color(const vec<3>& orig, const vec<3>& dir, vector<Sphere> spheres){
    double ans = numeric_limits<double>::max();
    // cout << ans << " ? ";
    double ans1 = ans;
    vec<3> res = bg_col;
    vec<3> N, point;
    for(int i = 0; i < spheres.size(); i++){
        // cout << spheres[i].center <<" "<<spheres[i].radius<<'\n';
        if(spheres[i].chk_ray(orig, dir, ans1, N) && ans1 < ans){
            // puts("QWQ");
            point = spheres[i].center + N * spheres[i].radius;
            // cout << N * N <<" ";
            // vec<3> light_dir = (light)
            // cout << (point - spheres[i].center).norm()<<" "<<spheres[i].radius;
            // float light_intensity
            ans = ans1, res = spheres[i].material.diffuse_color * get_diffuse_light(N, point) * spheres[i].material.albedo[0] 
                              + vec<3>{0, 255, 255, 255} * get_specular_light(N, point, (vec<3>{0, 0, 0, 0} - spheres[i].center).normalize(), spheres[i].material.specular_exponent) * spheres[i].material.albedo[1];
        }
    }
    res[1] = floor(res[1]), res[2] = floor(res[2]), res[3] = floor(res[3]);
    // if(res != bg_col)
    //     cout << res <<"\n";
    //     system("pause");
    // }
    return res;
    // if(sphere.chk_ray(orig, dir, ans)){
    //     return sphere.color;
    // }
    // return bg_col; 
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
            Image[i + (j - 1) * width] = (get_color(Camera, Now, spheres));
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
    Material m1(vec<2>{0.6, 0.3}, vec<3>{0, 124, 100, 32}, 50);
    Material m2(vec<2>{0.9, 0.1}, vec<3>{0, 5, 114, 107}, 10);
    // vec<3> co1 = vec<3>{0, 124, 100, 0};
    lights.push_back(Light(vec<3>{0, -20, 20, 20}, 0.5));
    lights.push_back(Light(vec<3>{0, 30, 50, -25}, 0.9));
    lights.push_back(Light(vec<3>{0, 30, 20, 30},  0.7));
    // co1[1] = 76, co1[2] = 208, co1[3] = 212;
    // ct[1] = 0, ct[2] = 0, ct[3] = -10;
    vector<Sphere> spheres;
    // vector<Material> materials;
    // spheres.push_back(Sphere(ct, 6, co1));
    spheres.push_back(Sphere(vec<3>{0, -3, 0, -16}, 2, m1));
    spheres.push_back(Sphere(vec<3>{0, -1, -1.5, -12}, 2, m2));
    spheres.push_back(Sphere(vec<3>{0, 1.5, -0.5, -18}, 3, m2));
    spheres.push_back(Sphere(vec<3>{0, 7, 5, -18}, 4, m1));
    work(spheres);
}