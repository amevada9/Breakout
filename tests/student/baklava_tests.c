#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "body.h"
#include "forces.h"
#include "list.h"
#include "polygon.h"
#include "scene.h"
#include "test_util.h"
#include "vector.h"

// 1 Test for each force

const double CLOSE_ENOUGH = .01;

list_t *make_shape() {
    size_t initial = 4;
    list_t *shape = list_init(initial, free);
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

bool close_to(vector_t vec1, vector_t vec2) {
    if ((vec1.x - CLOSE_ENOUGH < vec2.x) && (vec1.x + CLOSE_ENOUGH > vec2.x)) {
        return ((vec1.y - CLOSE_ENOUGH < vec2.y) && (vec1.y + CLOSE_ENOUGH > vec2.y));
    }
    return false;
}

void test_gravity() {
    double mass1 = 10;
    double mass2 = 20;
    double grav = 300;
    double timestep = 1;
    int steps = 5;
    scene_t *scene = scene_init();
    body_t *body1 = body_init(make_shape(), mass1, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, body1);
    body_set_centroid(body1, (vector_t){0, 0});
    body_t *body2 = body_init(make_shape(), mass2, (rgb_color_t){0, 0, 0});
    body_set_centroid(body2, (vector_t){50, 70});
    scene_add_body(scene, body2);
    create_newtonian_gravity(scene, grav, body1, body2);
    vector_t *vec_list_1 = malloc(sizeof(vector_t) * steps);
    assert(vec_list_1 != NULL && "The required memory was not allocated");
    vector_t *vec_list_2 = malloc(sizeof(vector_t) * steps);
    assert(vec_list_2 != NULL && "The required memory was not allocated");
    // actual gravity values (calculated)
    vec_list_1[0] = (vector_t){.2346, .3299};
    vec_list_1[1] = (vector_t){.9459, 1.3242};
    vec_list_1[2] = (vector_t){2.1448, 3.0027};
    vec_list_1[3] = (vector_t){3.8624, 5.4074};
    vec_list_1[4] = (vector_t){6.1507, 8.6110};
    vec_list_2[0] = (vector_t){49.8822, 69.8350};
    vec_list_2[1] = (vector_t){49.5270, 69.3379};
    vec_list_2[2] = (vector_t){48.9276, 68.4986};
    vec_list_2[3] = (vector_t){48.0688, 67.2963};
    vec_list_2[4] = (vector_t){46.9246, 65.6945};
    for (int i = 0; i < steps; i++) {
        scene_tick(scene, timestep);
        vector_t vector_1 = vec_list_1[i];
        vector_t vector_2 = vec_list_2[i];
        assert(close_to(vector_1, body_get_centroid(body1)));
        assert(close_to(vector_2, body_get_centroid(body2)));
    }
    free(vec_list_1);
    free(vec_list_2);
    scene_free(scene);
}

void test_drag() {
    double mass1 = 10;
    double gamma = 0.5;
    double timestep = 1;
    int steps = 5;
    scene_t *scene = scene_init();
    body_t *body1 = body_init(make_shape(), mass1, (rgb_color_t){0, 0, 0});
    vector_t vel1 = (vector_t){100, 100};
    body_set_velocity(body1, vel1);
    scene_add_body(scene, body1);
    create_drag(scene, gamma, body1);
    vector_t *vec_list_1 = malloc(sizeof(vector_t) * steps);
    assert(vec_list_1 != NULL && "The required memory was not allocated");
    // Test vectors calculated by hand
    vec_list_1[0] = (vector_t){97.5, 97.5};
    vec_list_1[1] = (vector_t){190.125, 190.125};
    vec_list_1[2] = (vector_t){278.11875, 278.11875};
    vec_list_1[3] = (vector_t){361.7128125, 361.7128125};
    vec_list_1[4] = (vector_t){441.1271719, 441.1271719};
    // Compare test vectors to the program
    for (int i = 0; i < steps; i++) {
        scene_tick(scene, timestep);
        vector_t vector_1 = vec_list_1[i];
        assert(close_to(vector_1, body_get_centroid(body1)));
    }
    free(vec_list_1);
    scene_free(scene);
}

bool passed(double past, double current, double a) {
    return ((past < a && current > a) || (past > a && current < a));
}

void test_spring_period() {
    const double M = 2;    // mass
    const double K = 5000; // k constant
    const double A = 15;   // distance a
    const double DT = 1e-6;
    const int STEPS = 1000000;
    const double period = 2 * M_PI * sqrt(M / K);
    // creating scene variables
    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t){A, 0});
    scene_add_body(scene, mass);
    body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, anchor);
    create_spring(scene, K, mass, anchor);
    int counter = 0;
    double time_counter = 0.0;
    double past = body_get_centroid(mass).x;
    double current;
    // running code
    for (int i = 0; i < STEPS; i++) {
        scene_tick(scene, DT);
        current = body_get_centroid(mass).x;
        if (passed(past, current, A)) {
            counter++;
            if (counter == 2) {
                assert(close_to((vector_t){time_counter, 0.0}, (vector_t){period, 0.0}));
                scene_free(scene);
                return;
            }
        }
        past = current;
        time_counter += DT;
    }
}

int main(int argc, char *argv[]) {
    // Run all tests if there are no command-line arguments
    bool all_tests = argc == 1;
    // Read test name from file
    char testname[100];
    if (!all_tests) {
        read_testname(argv[1], testname, sizeof(testname));
    }

    DO_TEST(test_gravity)
    DO_TEST(test_drag)
    DO_TEST(test_spring_period)

    puts("student_tests PASS");
}
