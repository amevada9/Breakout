#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "forces.h"
#include "test_util.h"
#include "vector.h"
#include "color.h"
#include "list.h"

/*** UTILITY FUNCTIONS ***/

const rgb_color_t COLOR = {0, 0, 0};

double random_double(void) {
    return ((double) rand() / RAND_MAX);
}

double random_double_range(double min, double max) {
    return min + random_double() * max;
}

vector_t random_vector_rect(double x_min, double x_max, double y_min, double y_max) {
    return (vector_t) {random_double_range(x_min, x_max),
                       random_double_range(y_min, y_max)};
}

vector_t random_vector_circ(double magnitude_min, double magnitude_max) {
    double angle = random_double_range(0, 2 * M_PI);
    double magnitude = random_double_range(magnitude_min, magnitude_max);
    return vec_rotate(vec_multiply(magnitude, (vector_t) {1, 0}), angle);
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

/*** TEST FOR EQUILIBRIUM POINTS ***/

void test_equil_gravity() {
    const size_t MIN_PTS = 2;
    const size_t MAX_PTS = 10;
    const double RADIUS  = 100;
    const int    TEST_MASS = 1;
    const int    ANCH_MASS = 100000;
    const double G       = 1000;
    const double DT      = 1e-5;
    const size_t STEPS   = 10000;

    for (size_t num_pts = MIN_PTS; num_pts <= MAX_PTS; num_pts++) {
        scene_t *scene = scene_init();

        // Add test particle.
        body_t *test_particle = body_init(make_shape(), TEST_MASS, COLOR);
        body_set_centroid(test_particle, VEC_ZERO);
        scene_add_body(scene, test_particle);

        // Add field sources.
        for (size_t i = 0; i < num_pts; i++) {
            body_t *anchor = body_init(make_shape(), ANCH_MASS, COLOR);
            vector_t position = vec_multiply(RADIUS, vec_rotate((vector_t) {1, 0},
                                             i * 2 * M_PI / num_pts));
            body_set_centroid(anchor, position);
            scene_add_body(scene, anchor);

            create_newtonian_gravity(scene, G, test_particle, anchor);
        }

        // Run simulation.
        for (size_t t = 0; t < STEPS; t++) {
            assert(vec_isclose(body_get_centroid(test_particle), VEC_ZERO));
            scene_tick(scene, DT);
        }

        scene_free(scene);
    }
}

void test_equil_spring() {
    const size_t MIN_PTS = 2;
    const size_t MAX_PTS = 10;
    const double RADIUS  = 100;
    const int    TEST_MASS = 1;
    const int    ANCH_MASS = 100000;
    const double K       = 1000;
    const double DT      = 1e-5;
    const size_t STEPS   = 10000;

    for (size_t num_pts = MIN_PTS; num_pts <= MAX_PTS; num_pts++) {
        scene_t *scene = scene_init();

        // Add test particle.
        body_t *test_particle = body_init(make_shape(), TEST_MASS, COLOR);
        body_set_centroid(test_particle, VEC_ZERO);
        scene_add_body(scene, test_particle);

        // Add field sources.
        for (size_t i = 0; i < num_pts; i++) {
            body_t *anchor = body_init(make_shape(), ANCH_MASS, COLOR);
            vector_t position = vec_multiply(RADIUS, vec_rotate((vector_t) {1, 0},
                                             i * 2 * M_PI / num_pts));
            body_set_centroid(anchor, position);
            scene_add_body(scene, anchor);

            create_spring(scene, K, test_particle, anchor);
        }

        // Run simulation.
        for (size_t t = 0; t < STEPS; t++) {
            assert(vec_isclose(body_get_centroid(test_particle), VEC_ZERO));
            scene_tick(scene, DT);
        }

        scene_free(scene);
    }
}

void test_equil() {
    test_equil_gravity();
    test_equil_spring();
}

/*** TEST FOR MOMENTUM CONSERVATION ***/

vector_t calculate_momentum(scene_t *scene) {
    vector_t p = VEC_ZERO;
    for (size_t i = 0; i < scene_bodies(scene); i++) {
        body_t *body = scene_get_body(scene, i);
        p = vec_add(p, vec_multiply(body_get_mass(body), body_get_velocity(body)));
    }
    return p;
}

double calculate_angular_momentum(scene_t *scene) {
    double l = 0;
    for (size_t i = 0; i < scene_bodies(scene); i++) {
        body_t *body = scene_get_body(scene, i);
        l += body_get_mass(body) * vec_cross(body_get_centroid(body),
                                             body_get_velocity(body));
    }
    return l;
}

void test_momentum_gravity() {
    const size_t NUM_PTS  = 50;
    const int    TEST_MASS = 10;
    const double INIT_VEL  = 100;
    const double INIT_POS  = 1000;
    const double G       = 1000;
    const double DT      = 1e-6;
    const size_t STEPS   = 10000;
    const double THRES   = 1e-3;

    scene_t *scene = scene_init();

    // Add bodies.
    for (size_t i = 0; i < NUM_PTS; i++) {
        body_t *body = body_init(make_shape(), TEST_MASS, COLOR);
        body_set_centroid(body, random_vector_rect(-INIT_POS, INIT_POS, -INIT_POS, INIT_POS));
        body_set_velocity(body, random_vector_circ(0, INIT_VEL));
        scene_add_body(scene, body);
    }

    // Add forces.
    for (size_t i = 0; i < NUM_PTS; i++) {
        for (size_t j = i + 1; j < NUM_PTS; j++) {
            create_newtonian_gravity(scene, G, scene_get_body(scene, i),
                                               scene_get_body(scene, j));
        }
    }

    vector_t p_init = calculate_momentum(scene);
    double   l_init = calculate_angular_momentum(scene);

    for (size_t t = 0; t < STEPS; t++) {
        assert(vec_within(THRES, calculate_momentum(scene), p_init));
        assert(within(THRES, calculate_angular_momentum(scene), l_init));
        scene_tick(scene, DT);
    }

    scene_free(scene);
}

void test_momentum() {
    test_momentum_gravity();
}

/*** TEST FOR DAMPING MONOTONICITY ***/

void test_damping_monotone() {
    const double GAMMA_MIN = 0.1;
    const double GAMMA_MAX = 100;
    const int    NUM_PTS   = 100;
    const int    TEST_MASS = 10;
    const double INIT_VEL  = 1000;
    const double DT      = 1e-6;
    const size_t STEPS   = 10000;

    scene_t *scene = scene_init();

    // Add bodies.
    for (size_t i = 0; i < NUM_PTS; i++) {
        body_t *body = body_init(make_shape(), TEST_MASS, COLOR);
        body_set_centroid(body, VEC_ZERO);
        body_set_centroid(body, (vector_t) {0, INIT_VEL});
        scene_add_body(scene, body);
        create_drag(scene, GAMMA_MIN + i * (GAMMA_MAX - GAMMA_MIN) / (NUM_PTS - 1), body);
    }

    for (size_t t = 0; t < STEPS; t++) {
        for (size_t i = 0; i < NUM_PTS - 1; i++) {
            assert(body_get_velocity(scene_get_body(scene, i)).y >=
                   body_get_velocity(scene_get_body(scene, i + 1)).y);
        }
        scene_tick(scene, DT);
    }

    scene_free(scene);
}

/*** TEST ZERO DAMPING COEFFICIENT USING SPRING FORCE ***/

void test_damping_spring_zero() {
    const int    TEST_MASS = 10;
    const double INIT_POS  = 1000;
    const double K         = 100;
    const double DT      = 1e-6;
    const size_t STEPS   = 100000;

    scene_t *scene = scene_init();

    body_t *anchor = body_init(make_shape(), INFINITY, COLOR);
    scene_add_body(scene, anchor);

    body_t *body1 = body_init(make_shape(), TEST_MASS, COLOR);
    body_t *body2 = body_init(make_shape(), TEST_MASS, COLOR);
    body_set_centroid(body1, (vector_t) {0, INIT_POS});
    body_set_centroid(body2, (vector_t) {0, INIT_POS});
    scene_add_body(scene, body1);
    scene_add_body(scene, body2);

    create_spring(scene, K, anchor, body1);
    create_spring(scene, K, anchor, body2);
    create_drag(scene, 0, body1);

    for (size_t t = 0; t < STEPS; t++) {
        assert(vec_isclose(body_get_centroid(body1), body_get_centroid(body2)));
        assert(vec_isclose(body_get_velocity(body1), body_get_velocity(body2)));
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

    DO_TEST(test_equil)
    DO_TEST(test_momentum)
    DO_TEST(test_damping_monotone)
    DO_TEST(test_damping_spring_zero)

    puts("student_tests PASS");
}
