#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "list.h"
#include "forces.h"
#include "test_util.h"
#include "vector.h"

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

/*
 * Tests that gravitational force between two objects, with the orbiting
 * object moving at velocity sqrt(G * (M1 + M2) / R) yields an orbit with
 * period 2 * pi * sqrt(R^3 / (G * (M1 + M2)))
 */
void test_gravitational_orbit() {
    const double M1 = 99999.9, M2 = 0.1;
    const double G = pow(2 * M_PI, 2);
    const double DT = 1e-5;
    const double R = 1e3;
    const int STEPS = 2500000;

    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass2, (vector_t){R, 0});
    body_set_velocity(mass2, (vector_t){0, sqrt(G * (M1 + M2) / R)});
    scene_add_body(scene, mass2);
    create_newtonian_gravity(scene, G, mass1, mass2);
    for (int i = 0; i < STEPS; i++) {
        scene_tick(scene, DT);
    }
    vector_t final_centroid = body_get_centroid(mass2);
    assert(within(1e-3, 0, final_centroid.x));
    assert(within(1e-3, R, final_centroid.y));
    scene_free(scene);
}

/*
 * Tests that a body decelerating due to drag has velocity given by
 * V0 * exp(-gamma * t / m)
 */
void test_drag_deceleration() {
    const double M = 1.7;
    const double V0 = 2.88;
    const double GAMMA = 1.25;
    const double DT = 1e-6;
    const int STEPS = 1000000;

    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t){0, 0});
    body_set_velocity(mass, (vector_t){V0, 0});
    scene_add_body(scene, mass);
    create_drag(scene, GAMMA, mass);
    for (int i = 0; i < STEPS; i++) {
        vector_t real = body_get_velocity(mass);
        vector_t expected = (vector_t) {V0 * exp(-GAMMA * i * DT / M), 0};
        assert(within(1e-5, real.x, expected.x));
        assert(real.y == expected.y);
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

/*
 * Tests that a mass on a spring with drag behaves like a damped harmonic
 * oscillator (underdamped case)
 */
void test_damped_oscillator() {
    const double M = 2;
    const double K = 3;
    const double GAMMA = 0.5;
    const double OMEGA = sqrt(4 * M * K - pow(GAMMA, 2)) / (2 * M);
    const double A = 2;
    const double B = GAMMA * A / (2 * M * OMEGA);
    const double AMP = sqrt(pow(A, 2) + pow(B, 2));
    const double PHI = atan(A / B);
    const double DT = 1e-6;
    const int STEPS = 1000000;

    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t){A, 0});
    scene_add_body(scene, mass);
    body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, anchor);
    create_spring(scene, K, mass, anchor);
    create_drag(scene, GAMMA, mass);
    for (int i = 0; i < STEPS; i++) {
        vector_t real = body_get_centroid(mass);
        vector_t expected = (vector_t) {
            AMP * exp(- GAMMA * i * DT / (2 * M)) * sin(OMEGA * i * DT + PHI),
            0
        };
        assert(within(1e-5, real.x, expected.x));
        assert(real.y == expected.y);
        assert(vec_equal(body_get_centroid(anchor), VEC_ZERO));
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

    DO_TEST(test_gravitational_orbit);
    DO_TEST(test_drag_deceleration);
    DO_TEST(test_damped_oscillator);

    puts("student_tests PASS");
}
