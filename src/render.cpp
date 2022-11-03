#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

int depth = 5;
const int numRays = 8; 
//Implementing the recursive ray-tracer 
//Source: Chapter 4.8, Fundamentals of Computer Graphics, Fourth Edition.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);
        Ray reflection = computeReflectionRay(ray, hitInfo);
        drawRay(ray, Lo);

        if (features.enableRecursive) {
            if (rayDepth >= depth || hitInfo.material.ks == glm::vec3{ 0.0 }) {
                return Lo;
            }
           
            if (features.extra.enableGlossyReflection) {
            
                return Lo + hitInfo.material.ks * glossyReflection(scene, bvh, reflection, 
                    features, rayDepth, hitInfo,
                    numRays);
            }
            //1/shininess
            //empty vector color 
            // for loop 
            //orthornomal basis 
            // epsilon, epsilon dash 
            //divide by number of sample 

            return Lo += hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
        }
     


        // Visual Debug: Draw a ray with a color which is the returned value from computeLightContribution
  
        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}


glm::vec3 glossyReflection(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth, HitInfo hitInfo, int numRays)
{
    float a = 1.0f / hitInfo.material.shininess;
    glm::vec3 color { 0.0f };
    glm::vec3 w; 
    //creating orthornomal basis
    if (glm::abs(hitInfo.normal.x) > 0.99) {
        w = glm::vec3(0, 1, 0);
    } else {
        w = glm::vec3(1, 0, 0); 
    }
    glm::vec3 u_vector = glm::normalize(glm::cross(hitInfo.normal, w));
    glm::vec3 v_vector = glm::normalize(glm::cross(hitInfo.normal, u_vector));

    for (int i = 0; i < numRays; i++) {
        for (int j = 0; j < numRays; j++) {
            float epsilon = std::rand() / RAND_MAX;
            float epsilon_dash = std::rand() / RAND_MAX;

            float u = -a / 2 + epsilon * a;
            float v = -a / 2 + epsilon_dash * a;

            glm::vec3 r_dash = ray.direction + u_vector * u + v_vector * v;

            glm::vec3 r_dash_origin = ray.origin + ray.t * ray.direction;

            float weight = glm::pow(glm::dot(glm::normalize(ray.direction), glm::normalize(r_dash)), hitInfo.material.shininess);

            color += getFinalColor(scene, bvh ,Ray {r_dash_origin, r_dash}, features, rayDepth + 1);
        }
    }
    color /= (numRays * numRays);
    return color;
}

void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    // Enable multi threading in Release mode
#ifdef NDEBUG
#pragma omp parallel for schedule(guided)
#endif
    for (int y = 0; y < windowResolution.y; y++) {
        for (int x = 0; x != windowResolution.x; x++) {
            // NOTE: (-1, -1) at the bottom left of the screen, (+1, +1) at the top right of the screen.
            const glm::vec2 normalizedPixelPos {
                float(x) / float(windowResolution.x) * 2.0f - 1.0f,
                float(y) / float(windowResolution.y) * 2.0f - 1.0f
            };
            const Ray cameraRay = camera.generateRay(normalizedPixelPos);
            screen.setPixel(x, y, getFinalColor(scene, bvh, cameraRay, features));
        }
    }
}