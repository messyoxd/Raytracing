#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "vector.hpp"

class Sphere{

    private:
        Vec3f center;
        float radius;
    public:
        Sphere(const Vec3f &c, const float &r): center(c), radius(r){}

        bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
            Vec3f L = center - orig;
            float tca = L*dir;
            float d2 = sqrt(L*L - tca*tca);
            if (d2 > radius*radius) return false;
            float thc = sqrtf(radius*radius - d2);
            t0       = tca - thc;
            float t1 = tca + thc;
            if (t0 < 0) t0 = t1;
            if (t0 < 0) return false;
            return true;
        }
};

Vec3f castRay(const Vec3f &orig, const Vec3f &dir, const Sphere &sphere){
    float sphere_dist = std::numeric_limits<float>::max();
    if(!sphere.ray_intersect(orig, dir, sphere_dist)){
        return Vec3f(0.8,0.7,0.8);
    }
    return Vec3f(0.4,0.4,0.3);
}

void render(const Sphere &sphere) {
    const int width    = 1024;
    const int height   = 768;
    const float fov = 1;
    std::vector<Vec3f> framebuffer(width*height);
    // obliquous case
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float u = (2*(j - 0.5)/(float) width-1)*(tan(fov/2.)*width/(float)height);
            float v = (2*(i-0.5)/(float) height-1)*(tan(fov/2.));
            Vec3f dir = Vec3f(u,v,-1).normalize();
            framebuffer[i+j*width] = castRay(Vec3f(0,0,0), dir, sphere);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Sphere sphere(Vec3f(0, 1, -8), 2);
    render(sphere);
    return 0;
}
