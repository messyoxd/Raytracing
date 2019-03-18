#include <limits>
#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include "vector.hpp"

class Light {

    public:
        Vec3f position;
        float intensity;
        Light(const Vec3f &p, const float &i): position(p), intensity(i){}

};

class Material {

    public:
        Vec3f diffuse_color;
        Material(const Vec3f &color): diffuse_color(color){}
        Material(): diffuse_color(){}

};

class Sphere{

    public:
        Vec3f center;
        float radius;
        Material material;

        Sphere(const Vec3f &c, const float &r, const Material &m): center(c), radius(r), material(m){}

        bool ray_intersect(const Vec3f &orig, const Vec3f &dir, float &t0) const {
            Vec3f L = center - orig;
            float tca = L*dir;
            if(tca < 0) return false;
            float d2 = L*L - tca*tca;
            if (d2 > radius*radius) return false;
            float thc = sqrt(radius*radius - d2);
            // t0 will be the distance between the center of the sphere
            // and the first point where the ray intersects
            t0       = tca - thc;
            float t1 = tca + thc;
            if (t0 < 0.0) t0 = t1;
            if (t0 < 0.0) return false;
            return true;
        }
};

/*
 * checks if the ray intersects with a object and determines what object
 * will interfere in the pixel's color (the closest to the camera)
*/
bool scene_intersect(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres,
                    Vec3f &hit, Vec3f &N, Material &material)
{
    float spheres_dist = std::numeric_limits<float>::max();
    for(size_t i=0; i < spheres.size(); i++){
        float dist_i;
        // dist_i will receive a value in ray_intersect
        if(spheres[i].ray_intersect(orig, dir, dist_i) && dist_i < spheres_dist){
            spheres_dist = dist_i;
            // hit is a vector that goes from the origin to the point where the ray intersects
            // the sphere
            hit = orig + dir*dist_i;
            // N is a directional unitary vector
            N = (hit - spheres[i].center).normalize();
            material = spheres[i].material;
        }
    }
    return spheres_dist<1000;
}

/*
 * Casts a ray in a pixel determined by dir
 */
Vec3f castRay(const Vec3f &orig, const Vec3f &dir, const std::vector<Sphere> &spheres,
              const std::vector<Light> &lights)
{
    Vec3f point, N;
    Material material;
    if(!scene_intersect(orig, dir, spheres, point, N, material)){
        return Vec3f(0.4,0.5,0.3);
    }
    float diffuse_light_intensity = 0.0;
    for (size_t i=0; i < lights.size(); i++){

        Vec3f light_dir = (lights[i].position - point).normalize();
        diffuse_light_intensity += (light_dir*N)*lights[i].intensity;

    }
    return material.diffuse_color*diffuse_light_intensity;
}

/*
 *  renders a scene
 */
void render(const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int width    = 1024;
    const int height   = 768;
    // const float fov = 3.14159265358979323846/2;
    std::vector<Vec3f> framebuffer(width*height);
    // obliquous case
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float u = (2*(i + 0.5)/(float) width-1);
            float v = -(2*(j+ 0.5)/(float) height-1);
            Vec3f dir = Vec3f(u,v,-1).normalize();
            framebuffer[i+j*width] = castRay(Vec3f(0,0,0), dir, spheres, lights);
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
    Material metal(Vec3f(0.2,0.3, 0.6));
    Material gold(Vec3f(0.99, 0.8, 0));
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3,0, -16),2,gold));
    spheres.push_back(Sphere (Vec3f(-1.0,-1.5,-18),2,metal));
    spheres.push_back(Sphere (Vec3f(3.0,-0.5,-15),2,metal));
    spheres.push_back(Sphere (Vec3f(7,5,-18),2,gold));
    std::vector<Light> lights;
    lights.push_back(Light (Vec3f(20,20,20),1.5));
    render(spheres, lights);
    return 0;
}
