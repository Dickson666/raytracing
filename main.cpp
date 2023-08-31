#include<cstdio>
#include<vector>
#include<cmath>
#include<iostream>
#include<fstream>
#include<limits>
#include"geometry.h"

// #define double float

using namespace std;

vec<3> bg_col;
struct Sphere{
    vec<3> center;
    double radius;
    vec<3> color;

    Sphere(const vec<3>& c, const double& r, const vec<3>& co) : center(c), radius(r), color(co) {}

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

double get_light(const vec<3>& N, const vec<3>& point){
    double diffuse_light_intensity = 0;
    for(int i = 0; i < lights.size(); i++){
        vec<3> light_dir = (lights[i].position - point).normalize();
        // cout << light_dir * N <<"\n";
        diffuse_light_intensity += lights[i].intensity * max(0.0, light_dir * N);
    }
    return diffuse_light_intensity;
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
            ans = ans1, res = spheres[i].color * get_light(N, point);
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
    int height = 720;
    int width = 1280;
    double fov = 2;

    vec<3> Camera;

    Camera[1] = 0, Camera[2] = 0, Camera[3] = 0;

    vector<vec<3> >Image(width * height);

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
    ofs << "P6\n" << width <<" "<< height<<"\n255\n";
    for(int i = 1; i <= width * height; i++){
        for(int j = 1; j <= 3; j++){
            // if(j != 3 && (abs(Image[i][j]) > 1e-5  && abs(Image[i][j] - 255) > 1e-5))
            //     cout << i <<" "<< j <<" "<< Image[i][j]<<" "<<Image[i][j] - 255<<" ",puts("QWQ");
            ofs << char(floor(Image[i][j]));
        }
        // ofs <<" , ";
        // if(i % width == 0)
        //     puts("");
    }
    ofs.close();
}
signed main(){
    bg_col[1] = bg_col[2] = bg_col[3] = 255;
    vec<3> ct = vec<3>{0, 0, 0, -10};
    vec<3> co1 = vec<3>{0, 0, 200, 0};
    lights.push_back(Light(vec<3>{0, 10, 10, 0}, 1));
    // co1[1] = 76, co1[2] = 208, co1[3] = 212;
    // ct[1] = 0, ct[2] = 0, ct[3] = -10;
    vector<Sphere> spheres;
    spheres.push_back(Sphere(ct, 6, co1));
    spheres.push_back(Sphere(vec<3>{0, 3, 5, -10}, 4, co1));
    work(spheres);
}
