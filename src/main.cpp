#include <cstdlib>
#include <cstdio>
#include <fstream>
#include <cstdint>
#include <bit>
#include <cmath>
#include <limits>
#include <random>
 
float f_inv(float n) {
  const float x2 = n * 0.5f;
  const float threehalfs = 1.5f;
  union { float f; uint32_t i; }
  conv = { .f = n };
  conv.i = 0x5f3759df - (conv.i >> 1);
  conv.f *= threehalfs - (x2 * conv.f * conv.f);
  return conv.f;
}

float rand_f() {
  static std::uniform_real_distribution<float> distribution(0.0, 1.0);
  static std::mt19937 generator;
  return distribution(generator);
}

float rand_f(float min, float max) { return min + (max - min) * rand_f(); }

struct Vec3t {
  float x, y, z;
  Vec3t() : x(0), y(0), z(0) { };
  Vec3t(float _x, float _y, float _z) : x(_x), y(_y), z(_z) { };
  Vec3t(float _e) : x(_e), y(_e), z(_e) { };
  // Multiply vector by another vector.
  Vec3t operator * (const Vec3t& v) const { return Vec3t(x * v.x, y * v.y, z * v.z); }
  // Add two vectors together.
  Vec3t operator + (const Vec3t& v) const { return Vec3t(x + v.x, y + v.y, z + v.z); }
  // Subtract a vector from another vector.
  Vec3t operator - (const Vec3t& v) const { return Vec3t(x - v.x, y - v.y, z - v.z); }
  // Multiply vector by a scalar.
  Vec3t operator * (float scalar) const { return Vec3t(x * scalar, y * scalar, z * scalar); }
  Vec3t operator *= (float scalar) const { return Vec3t(x + (x * scalar), y + (y * scalar), z + (z * scalar)); }
  // Get cross-product of two vectors.
  Vec3t operator % (const Vec3t& v) const { return Vec3t(y * v.z - z * v.y, z * v.x - x * v.z, x * v.y - y * v.x); }
  // Normalize vector.
  // Vec3t& norm() { return *this = *this * f_inv(x * x + y * y + z * z); }
  // Return dot-product of two matrices.
  float operator & (const Vec3t& v) const { return x * v.x + y * v.y + z * v.z; }
  float mag() const { return std::sqrt(x * x + y * y + z * z); }
  float mag2() const { return x * x + y * y + z * z; }

  static void printv(const Vec3t& v) {
    printf("%d %d %d\n", static_cast<int>(v.x), static_cast<int>(v.y), static_cast<int>(v.z));
  }
  static Vec3t reflect(const Vec3t& n, const Vec3t& v) {
    return v - (n * (2 * (v & n) / (n & n)));
  }
  static Vec3t random() { return Vec3t(rand_f(), rand_f(), rand_f()); }
  static Vec3t random(float min, float max) {
    return Vec3t(rand_f(min, max), rand_f(min, max), rand_f(min, max))); }
  static Vec3t random_us() {
    while(1) {
      Vec3t p = random(-1, 1);
      if(p.mag2() >= 1) continue;
      return p;
    }
  }
};
typedef struct Vec3t v3f;
typedef struct Vec3t p3f;
typedef struct Vec3t rgb;

struct Vec3t norm(const struct Vec3t& v) { 
  float len = f_inv(v.x * v.x + v.y * v.y + v.z * v.z);
  return Vec3t(v.x * len, v.y * len, v.z * len); }

struct Cord {
  v3f o, d;
  Cord() : o(), d() { };
  Cord(v3f _o, v3f _d) : o(_o), d(_d) { }; 
  p3f at(float t) const { return o + (d * t); }
};
typedef struct Cord ray;

struct Light {
  p3f o;
  float intensity;
  Light() : o(), intensity() { };
  Light(v3f _o, float _int) : o(_o), intensity(_int) { };
};

struct Material {
  float rough_index, reflect_index, trans_index;
  Material() : rough_index(), reflect_index(), trans_index() { };
  Material(float _ro, float _re, float _tr) : rough_index(_ro), reflect_index(_re), trans_index(_tr) { };
}; typedef struct Material material;

struct Plane {
  p3f p1, p2, p3, p4;
  rgb color;
  Plane()
   : p1(), p2(), p3() { };
  Plane(p3f _p1, p3f _p2, p3f _p3)
   : p1(_p1), p2(_p2), p3(_p3) {
     p4 = p2 + p3;
  };
  v3f get_normal() {
    v3f v = p2 - p1;
    v3f u = p3 - p1;
    return u % v;
  }
  v3f get_center() {
    return (p2 * 0.5) + (p3 * 0.5);
  }
};

struct Sphere {
  p3f center;
  float radius;
  material mat;
  Sphere()
   : center(), radius(), mat() { };
  Sphere(p3f _c, float _r)
   : center(_c), radius(_r), mat() { };
  Sphere(p3f _c, float _r, struct Material _m)
   : center(_c), radius(_r), mat(_m) { }
};

float intersect_sphere(const ray& r, struct Sphere s) {
  v3f oc = r.o - s.center;
  float a = r.d & r.d;
  float b = 2.0 * (oc & r.d);
  float c = (oc & oc) - s.radius * s.radius;
  float disc = b * b - 4 * a * c;
  if(disc < 0) return -1.0;
  return (-b - sqrt(disc)) * (1 / a) * 0.5;
}

float intersect_plane(const ray& r, struct Plane p) {
   float denom = (p.get_normal() & r.d);
   if(fabs(denom) > 0.0001) {
     v3f diff = p.get_center() - r.o;
     float t = (diff & p.get_normal()) / denom;
     if(t > 0.0001)
      return t;
   }
   return -1.0;
}

static struct Plane planes[1] = {
  Plane(0, 0, 1)
};

static struct Sphere spheres[] = {
  Sphere(p3f(1, -0.5, -3), 1, material(0, 0.5, 0)),
  Sphere(p3f(0, -101.5, -3), 100),
  Sphere(p3f(-1, -0.5, -3), 1)
};

static struct Light light_sources[] = {
  Light(p3f(4, 10, -3), 0.7),
  Light(p3f(-4, 10, -3), 0.7)
};

float photon_cast(p3f at, int depth = 1) {
  float max_lint = 0.2;
  for(int light = 0; light < sizeof(light_sources) / sizeof(Light); light++) {
    int is_blocked = 0;
    for(int object = 0; object < sizeof(spheres) / sizeof(Sphere); object++) {
      if(intersect_sphere(ray(at, light_sources[light].o - at), spheres[object]) > 0.0) { is_blocked++; break; } 
    }    
    if(!is_blocked) 
      max_lint = (light_sources[light].intensity > max_lint) ? light_sources[light].intensity : max_lint;
  }
  return max_lint;
}

rgb ray_color(const ray& r, int depth = 1) {
  if(depth <= 0) return rgb(0, 0, 0);
  float t, t_temp;
  int has_hit = 0;
  v3f n, v, u_d = norm(r.d);
  t = 0.5 * (u_d.y + 1.0);
  rgb rcol = (rgb(255, 255, 255) * (1.0 - t)) + (rgb(128, 155, 255) * t); 
  for(int o = 0; o < sizeof(spheres) / sizeof(Sphere); o++) {
    Sphere obj = spheres[o];
    t = intersect_sphere(r, obj);
    if(t <= 0.0) continue;
    if(t_temp > t || !has_hit) { 
      has_hit++;
      v = r.at(t) - spheres[o].center;
      n = norm(v);
      float lint = photon_cast(r.at(t - 0.00001));
      rcol = (((rgb(n.x, n.y, n.z) * 255) + rgb(255, 255, 255)) * (0.5 * lint));
    }
    t_temp = t;
  }
  return rcol;
}

int main() {
  // Image
  const float aspect_ratio = 1.0 / 1.0;
  const int image_width = 400;
  const int image_height = static_cast<int>(image_width / aspect_ratio);

  // Camera
  const float vp_h = 2.0;
  const float vp_w = aspect_ratio * vp_h;
  const float focal_length = 1.0;

  p3f camera = p3f(0);
  v3f horizontal = v3f(vp_w, 0, 0);
  v3f vertical = v3f(0, vp_h, 0);
  p3f lower_left_corner = camera - (horizontal * 0.5) - (vertical * 0.5) - v3f(0, 0, focal_length);

  printf("P3\n%d %d\n255\n", image_width, image_height);
  for(int j = image_height - 1; j >= 0; j--) {
    for(int i = 0; i < image_width; i++) {
      double u = double(i) / (image_width - 1);
      double  v = double(j) / (image_height - 1);
      ray r(camera, lower_left_corner + (horizontal * u) + (vertical * v) - camera);
      rgb color = ray_color(r);

      Vec3t::printv(color);
    }
  }
}
