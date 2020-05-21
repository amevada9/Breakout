#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "forces.h"
#include "test_util.h"


const double ANGLE_BETWEEN_POINTS = 0.05;
const double EPSILON = 0.001;
const double ANGLE_OF_CIRCLE = 360;
const double G = 1e3;
const double DT = 1e-6;
const int STEPS = 100;
const double RADIUS = 4;

const double MASS_PLANET = 1000;
const double MASS_SATELLITE = 0.01;
const rgb_color_t COLOR_SATELLITE = {1, 0, 1};
const rgb_color_t COLOR_PLANET = {0, 1, 0};
const vector_t CENTROID_PLANET = {50.0, 50.0};
double ORBITAL_RADIUS;

list_t *circle_sector(size_t rad, double min_rad, double max_rad){
    vector_t center = VEC_ZERO;
    list_t *points = list_init(1, free);
    double theta = min_rad;
    if (fabs(min_rad - 0) > EPSILON || fabs(max_rad - 2 * M_PI) > EPSILON){
        vector_t *c = malloc(sizeof(vector_t));
        *c = center;
        list_add(points, c);
    }
    while(theta < max_rad){
        vector_t *vertex_add = malloc(sizeof(vector_t));
        vertex_add->x = cos(theta)*rad;
        vertex_add->y = sin(theta)*rad;
        vector_t *point = malloc(sizeof(vector_t));
        *point = vec_add(center, *vertex_add);
        list_add(points, point);
        theta += ANGLE_BETWEEN_POINTS;
        free(vertex_add);
    }
    return points;
}


list_t *make_shape() {
    list_t *shape = list_init(4, free);
    vector_t *v = malloc(sizeof(*v));
    *v = (vector_t){-1, -1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){+1, -1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){+1, +1};
    list_add(shape, v);
    v = malloc(sizeof(*v));
    *v = (vector_t){-1, +1};
    list_add(shape, v);
    return shape;
}

void check_radius(scene_t *sc, body_t *s) {
  for (int i = 0; i < STEPS; i++) {
      double dx = body_get_centroid(s).x - CENTROID_PLANET.x;
      double dy = body_get_centroid(s).y - CENTROID_PLANET.y;
      double dist_sq = dx * dx + dy * dy;

      double distance = sqrt(dist_sq);
      assert(within(1e-4, distance / ORBITAL_RADIUS, 1));
      scene_tick(sc, DT);
  }
}

void test_gravity_orbit() {
    scene_t *scene = scene_init();
    body_t *satellite = body_init(make_shape(), MASS_SATELLITE, COLOR_SATELLITE);
    scene_add_body(scene, satellite);
    body_t *planet = body_init(make_shape(), MASS_PLANET, COLOR_PLANET);
    body_set_centroid(planet, CENTROID_PLANET);
    scene_add_body(scene, planet);
    create_newtonian_gravity(scene, G, satellite, planet);
    vector_t curr_s = body_get_centroid(satellite);
    double dx = (curr_s.x - CENTROID_PLANET.x) * (curr_s.x - CENTROID_PLANET.x);
    double dy = (curr_s.y - CENTROID_PLANET.y) * (curr_s.y - CENTROID_PLANET.y);
    ORBITAL_RADIUS = sqrt(dx + dy);

    check_radius(scene, satellite);
    scene_free(scene);
}


body_t *make_anchor(vector_t center){
    list_t *anchor_circle = circle_sector(RADIUS, 0, 2*M_PI);
    body_t *anchor = body_init(anchor_circle, INFINITY, (rgb_color_t) {1, 1, 1});
    body_set_centroid(anchor, center);
    return anchor;
}

void vertical_spring_position(){
    const double Mearth = 5.98e24;
    const double Mweight = 100.0;
    const double G = 6.673e-11;
    const double A = -75.0;
    const double POSearth = -6.3800000000000e6;
    double POSeq = -20.0;
    double K = -G * Mweight * Mearth / (pow(POSearth - POSeq, 2) * POSeq);
    double delta;
    scene_t *scene = scene_init();
    body_t *earth = body_init(circle_sector(2.0, 0, 2*M_PI), Mearth, (rgb_color_t) {0, 1, 0});
    body_set_centroid(earth, (vector_t) {0, POSearth});
    scene_add_body(scene, earth);
    body_t *anchor= body_init(circle_sector(2, 0, 2*M_PI), INFINITY, (rgb_color_t) {1, 1, 1});
    body_set_centroid(anchor, VEC_ZERO);
    scene_add_body(scene, anchor);
    body_t *weight = body_init(make_shape(), Mweight, (rgb_color_t) {0, 0, 0});
    body_set_centroid(weight, (vector_t) {0, POSeq + A});
    scene_add_body(scene, weight);

    create_spring(scene, K, weight, anchor);
    create_newtonian_gravity(scene, G, weight, earth);

    for(int i = 0; i < STEPS; i++){
        vector_t w = body_get_centroid(weight);
        vector_t v = (vector_t) {0, A  *  cos(sqrt(K / Mweight) * i * DT) + POSeq};
        // margin for weight position dependency
        if (i < 29400){
            delta = 1e-7;
        }
        else if (i < 93240){
            delta = 1e-6;
        }
        else if (i < 296200){
            delta = 1e-5;
        }
        else if (i < 975000){
            delta = 1e-4;
        }
        else{
            delta = 1e-3;
        }
        assert(within(delta, w.y, v.y));
        //printf("Earth Centroid: %f, %f\n", body_get_centroid(earth).x, body_get_centroid(earth).y);
        //printf("Expected: %f, %f\n", 0.0, POSearth);
        assert(vec_isclose(body_get_centroid(earth), (vector_t) {0, POSearth}));
        //printf("Anchor Centroid: %f, %f\n", 100 * body_get_centroid(anchor).x, 100 * body_get_centroid(anchor).y);
        assert(vec_equal(body_get_centroid(anchor), VEC_ZERO));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

double spring_potential(double K, body_t *body1, body_t *body2){
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return 0.5 * K * vec_dot(r, r);
}

double kinetic_energy(body_t *body){
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

void test_spring_energy_conservation(){
    const double M = 10.0;
    const double K = 2.0;
    const double A = 3.0;
    scene_t *scene = scene_init();
    body_t *anchor= body_init(circle_sector(2, 0, 2*M_PI), INFINITY, (rgb_color_t) {1, 1, 1});
    body_set_centroid(anchor, VEC_ZERO);
    scene_add_body(scene, anchor);
    body_t *weight = body_init(make_shape(), M, (rgb_color_t) {0, 0, 0});
    body_set_centroid(weight, (vector_t) {A, 0});
    scene_add_body(scene, weight);
    create_spring(scene, K, weight, anchor);
    double initial_energy = spring_potential(K, weight, anchor);
    for (int i = 0; i < STEPS; i++){
        double energy = spring_potential(K, weight, anchor) + kinetic_energy(weight);
        assert(within(1e-7, energy / initial_energy, 1));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

int main(int argc, char *argv[]) {
    bool all_tests = argc == 1;
    char testname[100];
    if (!all_tests) {
        read_testname(argv[1], testname, sizeof(testname));
    }

    DO_TEST(test_gravity_orbit);
    DO_TEST(test_spring_energy_conservation);
    DO_TEST(vertical_spring_position);

    puts("forces_test PASS");
}
