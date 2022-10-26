#include "render.h"
#include "intersect.h"
#include "light.h"
#include "screen.h"
#include <framework/trackball.h>
#ifdef NDEBUG
#include <omp.h>
#endif

//Implementing the recursive ray-tracer 
//Source: Chapter 4.8, Fundamentals of Computer Graphics, Fourth Edition.
glm::vec3 getFinalColor(const Scene& scene, const BvhInterface& bvh, Ray ray, const Features& features, int rayDepth)
{
    HitInfo hitInfo;
    if (bvh.intersect(ray, hitInfo, features)) {

        glm::vec3 Lo = computeLightContribution(scene, bvh, features, ray, hitInfo);

        if (features.enableRecursive) {
            Ray reflection = computeReflectionRay(ray, hitInfo);
            if (rayDepth < 10 && hitInfo.material.ks != glm::vec3{ 0.0f }) {
                Lo += hitInfo.material.ks * getFinalColor(scene, bvh, reflection, features, rayDepth + 1);
            }
        }
        //if normal interpolation flag enables, we draw normal at 3 vertices and the interpolated normal.    
        if (features.enableNormalInterp) {
            drawRay(Ray { hitInfo.vertices[0].position, hitInfo.vertices[0].normal, 2 }, glm::vec3 { 1.0f, 0.0f, 0.0f });
            drawRay(Ray { hitInfo.vertices[1].position, hitInfo.vertices[1].normal, 2 }, glm::vec3 { 1.0f, 0.0f, 0.0f });
            drawRay(Ray { hitInfo.vertices[2].position, hitInfo.vertices[2].normal, 2 }, glm::vec3 { 1.0f, 0.0f, 0.0f });
            drawRay(Ray { ray.origin + ray.t * ray.direction,
                      hitInfo.normal, 2 },
             glm::vec3 { 0.0f, 1.0f, 0.0f });
        }


        // Visual Debug: Draw a ray with a color which is the returned value from computeLightContribution
        drawRay(ray, Lo);
  
        // Set the color of the pixel to white if the ray hits.
        return Lo;
    } else {
        // Draw a red debug ray if the ray missed.
        drawRay(ray, glm::vec3(1.0f, 0.0f, 0.0f));
        // Set the color of the pixel to black if the ray misses.
        return glm::vec3(0.0f);
    }
}

glm::vec3 pixelAt(int i, int j, Screen& screen) {
    if (i < 0 || i >= screen.resolution().x || j < 0 || j >= screen.resolution().y ) {
        return glm::vec3 { 0.0f };
    }
    return screen.pixels()[screen.indexAt(i, j)];
}

glm::vec3 boxFilter(int i, int j, int filterSize, Screen& screen, float scale) {
    filterSize = std::max(1, filterSize);
    glm::vec3 result { 0.0f };
    for (int x = -filterSize; x < filterSize + 1; ++x) {
        for (int y = -filterSize; y < filterSize + 1; ++y) {
            result += pixelAt(i, j, screen);
        }
    }
    result *= scale / ((2.0f * filterSize + 1.0f) * (2.0f * filterSize + 1.0f));
    return result;
}



void renderRayTracing(const Scene& scene, const Trackball& camera, const BvhInterface& bvh, Screen& screen, const Features& features)
{
    glm::ivec2 windowResolution = screen.resolution();
    int width = windowResolution.x;
    int height = windowResolution.y;
    glm::vec3 inital { 0.0f };
    std::vector<glm::vec3> pixels(width*height);

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
            glm::vec3 color = getFinalColor(scene, bvh, cameraRay, features);
            screen.setPixel(x, y, color);
          
        }
    }



     if (features.extra.enableBloomEffect) {
        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                glm::vec3 color = screen.pixels()[screen.indexAt(x, y)];
                if (color.x > 0.8 && color.y > 0.8 > color.z > 0.8) {
                    pixels[screen.indexAt(x, y)] = color;
                } else {
                    pixels[screen.indexAt(x, y)] = glm::vec3 { 0.0f };
                }
            }
        }
        

        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                pixels[screen.indexAt(x, y)] = boxFilter(x, y, 1, screen, 2);
            }
        }

        for (int y = 0; y < windowResolution.y; y++) {
            for (int x = 0; x != windowResolution.x; x++) {
                if (pixels[screen.indexAt(x, y)] != glm::vec3 { 0.0f }) {
                    screen.setPixel(x, y, (screen.pixels()[screen.indexAt(x, y)] + pixels[screen.indexAt(x, y)]) / 2.0f);
                }
            }
        }
    }
}