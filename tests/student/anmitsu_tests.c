#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "forces.h"
#include "test_util.h"

list_t *make_default_shape() {
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

void test_gravity_stationary_center() {
    const double M1 = 1.61803399, M2 = 1.61803399;
    const double G = 1e3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_default_shape(), M1, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass1, (vector_t){-5, 10});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_default_shape(), M2, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass2, (vector_t){5, -10});
    scene_add_body(scene, mass2);
    create_newtonian_gravity(scene, G, mass1, mass2);
    for (int i = 0; i < STEPS; i++) {
        assert(vec_isclose(vec_add(body_get_centroid(mass1), body_get_centroid(mass2)), VEC_ZERO));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

void test_spring_stationary_center() {
    const double M1 = 3.14*2.72, M2 = 3.14*2.72;
    const double K = 3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_default_shape(), M1, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass1, (vector_t){-10, -10});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_default_shape(), M2, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass2, (vector_t){10, 10});
    scene_add_body(scene, mass2);
    create_spring(scene, K, mass1, mass2);
    for (int i = 0; i < STEPS; i++) {
        assert(vec_isclose(vec_add(body_get_centroid(mass1), body_get_centroid(mass2)), VEC_ZERO));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

double spring_potential(double K, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return 0.5 * K * vec_dot(r, r);
}

double kinetic_energy(body_t *body) {
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

void test_spring_energy_conservation() {
    const double M1 = 2.72, M2 = 3.14;
    const double K = 2;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_default_shape(), M1, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_default_shape(), M2, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass2, (vector_t){-20, 10});
    scene_add_body(scene, mass2);
    create_spring(scene, K, mass1, mass2);
    double initial_energy = spring_potential(K, mass1, mass2);
    for (int i = 0; i < STEPS; i++) {
        double energy = spring_potential(K, mass1, mass2) + kinetic_energy(mass1) +
                        kinetic_energy(mass2);
        assert(within(1e-4, energy / initial_energy, 1));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

int main(int argc, char *argv[]) {
    // Run all tests if there are no command-line arguments
    bool all_tests = argc == 1;
    // Read test name from file
    char testname[100];
    if (!all_tests) {
        read_testname(argv[1], testname, sizeof(testname));
    }

    DO_TEST(test_spring_energy_conservation);
    DO_TEST(test_spring_stationary_center);
    DO_TEST(test_gravity_stationary_center);

    puts("student_test PASS");
}
