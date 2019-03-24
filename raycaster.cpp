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
        Vec4f albedo;
        float refractive_index;
        float specular_coeficient;
        Material(const Vec3f &color, const Vec4f &a, const float &spec, const float &r): diffuse_color(color), albedo(a), specular_coeficient(spec), refractive_index(r){}
        Material(): diffuse_color(), albedo(1,0,0,0), specular_coeficient(), refractive_index(1){}

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
 *  refract with Snell's law
*/
Vec3f refract(const Vec3f &I, const Vec3f &N, const float &refractive_index){

    float cosi = - std::max(-1.f, std::min(1.f, I*N));
    float etai = 1, etat = refractive_index;
    Vec3f n = N;
    if(cosi < 0){
        cosi = -cosi;
        std::swap(etai,etat); n = -N;
    }
    float eta = etai/etat;
    float k = 1 - eta*eta*(1 - cosi*cosi);
    return k < 0 ? Vec3f(0,0,0) : I*eta + n*(eta * cosi - sqrtf(k));

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
              const std::vector<Light> &lights, size_t depth=0)
{
    Vec3f point, N;
    Material material;
    if(depth>4 || !scene_intersect(orig, dir, spheres, point, N, material)){
        return Vec3f(0,0,0);
    }
    Vec3f reflect_dir = reflect(dir, N);
    Vec3f refract_dir = refract(dir, N, material.refractive_index).normalize();
    Vec3f reflect_orig = reflect_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f refract_orig = refract_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
    Vec3f reflect_color = castRay(reflect_orig, reflect_dir, spheres, lights, depth+1);
    Vec3f refract_color = castRay(refract_orig, refract_dir, spheres, lights, depth+1);

    float diffuse_light_intensity = 0.0, specular_light_intensity = 0;
    for (size_t i=0; i < lights.size(); i++){

        Vec3f light_dir = (lights[i].position - point).normalize();
        // shadows
        float light_distance = (lights[i].position - point).norm();
        Vec3f shadow_orig = light_dir*N < 0 ? point - N*1e-3 : point + N*1e-3;
        Vec3f shadow_pt, shadow_N;
        Material tmpmaterial;
        if (scene_intersect(shadow_orig, light_dir, spheres, shadow_pt, shadow_N, tmpmaterial) && (shadow_pt-shadow_pt).norm() < light_distance)
            continue;
        //Lambert's reflection model
        diffuse_light_intensity += (light_dir*N)*lights[i].intensity;

        // Phong's reflection model
        // specular_light_intensity = (0, R.V)^specular_index, where V is the vector
        // from the pixel to the object's intersection point and R is the vector that
        // is created by the reflected ray from the light source by the object's surface
        specular_light_intensity += powf(std::max(0.f, -reflect(-light_dir,N)*dir),material.specular_coeficient)*lights[i].intensity;

    }
    return material.diffuse_color*diffuse_light_intensity * material.albedo[0] + Vec3f(1.,1.,1.)*specular_light_intensity * material.albedo[1] + reflect_color*material.albedo[2] + refract_color*material.albedo[3];
}

/*
 *  renders a scene
 */
void render(Vec3f &camera, const std::vector<Sphere> &spheres, const std::vector<Light> &lights) {
    const int width    = 1000;
    // const int width    = 100;
    const int height   = 1000;
    // const int height   = 100;
    const float fov = 3.14159265358979323846/2;
    // const float fov = 2.8;
    std::vector<Vec3f> framebuffer(width*height);

    Vec3f w = camera.normalize();

    Vec3f t = w.getLowerComponent();
    printf("vetor w %f %f %f \n",w.getX(),w.getY(),w.getZ())
    printf("vetor t %f %f %f \n",t.getX(),t.getY(),t.getZ());
    Vec3f u = (cross(t,w)).normalize();
    printf("vetor u %f %f %f \n",u.getX(),u.getY(),u.getZ());
    Vec3f v = cross(u,w);
    printf("vetor v %f %f %f \n",v.getX(),v.getY(),v.getZ());

    // orthogonal case
    for (size_t j = 0; j<height; j++) {
        for (size_t i = 0; i<width; i++) {
            float u_linha = -5+(10*(i + 0.5)/(float) width-1)*tan(fov/2.)*width/(float)height;
            // float u = (2*(i + 0.5)/(float) width-1);
            float v_linha = -5+(10*(j+ 0.5)/(float) height-1)*tan(fov/2.);
            // float v = -(2*(j+ 0.5)/(float) height-1);
            Vec3f dir = w*-1;
            Vec3f origin = (camera + u*u_linha + v*v_linha);
            framebuffer[i+j*width] = castRay(origin, dir, spheres, lights);
        }
    }

    // // obliquous case
    // for (size_t j = 0; j<height; j++) {
    //     for (size_t i = 0; i<width; i++) {
    //         float u_linha = -5+(10*(i + 0.5)/(float) width-1)*tan(fov/2.)*width/(float)height;
    //         // float u = (2*(i + 0.5)/(float) width-1);
    //         float v_linha = -5+(10*(j+ 0.5)/(float) height-1)*tan(fov/2.);
    //         // float v = -(2*(j+ 0.5)/(float) height-1);
    //         Vec3f dir = (w*-1+ u*u_linha + v*v_linha).normalize();
    //         // Vec3f dir = Vec3f(-1,u_linha,v_linha).normalize();
    //         Vec3f origin = camera;
    //         framebuffer[i+j*width] = castRay(origin, dir, spheres, lights);
    //     }
    // }

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

    Material trabalho(Vec3f(0.2, 0.3, 0.6), Vec4f(1,1,0,0), 500, 1.0);
    // this is where the pixel mail is
    Vec3f camera = Vec3f(5,5,0);
    std::vector<Sphere> spheres;

    spheres.push_back(Sphere(Vec3f(0,0,0),1,trabalho));

    std::vector<Light> lights;
    
    lights.push_back(Light (Vec3f(4,4,4),1.5));
    render(camera, spheres, lights);
    return 0;
}
