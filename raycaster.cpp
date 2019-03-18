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
        Vec2f albedo;
        float specular_gradient;
        Material(const Vec3f &color, const Vec2f &a, const float &spec): diffuse_color(color), albedo(a), specular_gradient(spec){}
        Material(): diffuse_color(), albedo(1,0), specular_gradient(){}

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
 *  calculates the reflection ray
*/
Vec3f reflect(const Vec3f &I, const Vec3f &N){
    return I - N*2.f*(I*N);
}

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
        return Vec3f(1.,1.,1.);
    }
    float diffuse_light_intensity = 0.0, specular_light_intensity = 0;
    for (size_t i=0; i < lights.size(); i++){

        Vec3f light_dir = (lights[i].position - point).normalize();

        diffuse_light_intensity += (light_dir*N)*lights[i].intensity;
        // specular_light_intensity = (0, R.V)^specular_index, where V is the vector
        // from the pixel to the object's intersection point and R is the vector that
        // is created by the reflected ray from the light source by the object's surface
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir,N)*dir),material.specular_gradient)*lights[i].intensity;

    }
    return material.diffuse_color*diffuse_light_intensity * material.albedo[0] + Vec3f(1.,1.,1.)*specular_light_intensity * material.albedo[1];
}

/*
 *  renders a scene
 */
void render(const Vec3f &camera, const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
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
            framebuffer[i+j*width] = castRay(camera, dir, spheres, lights);
        }
    }

    std::ofstream ofs; // save the framebuffer to file
    ofs.open("./out.ppm");
    ofs << "P6\n" << width << " " << height << "\n255\n";
    for (size_t i = 0; i < height*width; ++i) {
        Vec3f &c = framebuffer[i];
        float max = std::max(c[0], std::max(c[1],c[2]));
        if(max>1) c=c*(1./max);
        for (size_t j = 0; j<3; j++) {
            ofs << (char)(255 * std::max(0.f, std::min(1.f, framebuffer[i][j])));
        }
    }
    ofs.close();
}

int main() {
    Material metal(Vec3f(0.44,0.44, 0.44), Vec2f(0.6,0.3), 50.);
    Material gold(Vec3f(0.99, 0.8, 0), Vec2f(0.6,0.3), 50.);
    Material red_rubber(Vec3f(0.3, 0.1, 0.1),Vec2f(0.9,  0.1),   10.);
    // this is where the pixel mail is
    Vec3f camera = Vec3f(0,0,0);
    std::vector<Sphere> spheres;
    spheres.push_back(Sphere(Vec3f(-3,0, -16),2,gold));
    spheres.push_back(Sphere (Vec3f(-1.0,-1.5,-18),2,metal));
    spheres.push_back(Sphere (Vec3f(3.0,-0.5,-15),2,metal));
    spheres.push_back(Sphere (Vec3f(7,5,-18),2,gold));
    spheres.push_back(Sphere (Vec3f(0,10,-18),2,red_rubber));
    std::vector<Light> lights;
    lights.push_back(Light (Vec3f(-20,20,20),1.5));
    lights.push_back(Light (Vec3f(30,-50,-25),1.8));
    lights.push_back(Light (Vec3f(30,20,30),1.7));
    render(camera, spheres, lights);
    return 0;
}
