#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include "list.h"
#include "vector.h"
#include "forces.h"
#include "body.h"
#include "test_util.h"

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

double kinetic_energy(body_t *body) {
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

double spring_potential(body_t *body1, body_t *body2, double k) {
    // Equilibrium is when they are completely overlapping
    // Assuming 1D case (x)
    vector_t p1 = body_get_centroid(body1);
    vector_t p2 = body_get_centroid(body2);
    vector_t diff = vec_subtract(p2, p1);
    return 0.5 * k * pow(diff.x,2);
}

double gravity_potential(double G, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return -G * body_get_mass(body1) * body_get_mass(body2) / sqrt(vec_dot(r, r));
}

void test_spring_conservation() { // Will only work if spring collision isn't a thing
    const double M = 10;
    const double K = 2;
    const double A = 3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t){A, 0});
    scene_add_body(scene, mass);
    body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, anchor);
    create_spring(scene, K, mass, anchor);

    double initial_energy = spring_potential(mass, anchor, K);
    for (int i = 0; i < STEPS; i++) {
        // Want to assert that the KE + PE = E
        double potential = spring_potential(mass, anchor, K);
        double kinetic = kinetic_energy(mass);
        assert(within(1e-4, (potential + kinetic) / initial_energy, 1));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

void test_fast_gravity_conservation(){ // Tests accuracy of the gravity
    const double M1 = 100, M2 = 50000;
    const double G = 1e6;
    const double DIST = 1e5;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass2, (vector_t){DIST, 0});
    scene_add_body(scene, mass2);
    create_newtonian_gravity(scene, G, mass1, mass2);
    double initial_energy = gravity_potential(G, mass1, mass2);
    for (int i = 0; i < STEPS; i++) {
        assert(body_get_centroid(mass1).x < body_get_centroid(mass2).x); // ???
        double energy = gravity_potential(G, mass1, mass2) + kinetic_energy(mass1) +
                        kinetic_energy(mass2);
        assert(within(1e-8, energy / initial_energy, 1));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

void test_drag_slows(){
    const double M1 = 10, M2 = 10, M3 = 1000;
    const double G = 10;
    const double DIST = 10.0;
    const double DT = 1e-6;
    const int STEPS = 10000;
    const double GAMMA = 1;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_shape(), M1, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_shape(), M2, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, mass2);
    body_t *mass3 = body_init(make_shape(), M3, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass3, (vector_t){DIST, 0});
    scene_add_body(scene, mass3);

    create_newtonian_gravity(scene, G, mass1, mass3);
    create_newtonian_gravity(scene, G, mass2, mass3);
    create_drag(scene, GAMMA, mass2);

    for (int i = 0; i < STEPS; i++) {
        double energy1 = fabs(gravity_potential(G, mass1, mass3)) + kinetic_energy(mass1);
        double energy2 = fabs(gravity_potential(G, mass2, mass3)) + kinetic_energy(mass2);
        if (i > 100) {
            assert(energy1 > energy2);
        }
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

void test_spring_center(){
    const double M = 10;
    const double K = 2;
    const double A = 4;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass1 = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass1, (vector_t){-A/2, 0});
    scene_add_body(scene, mass1);
    body_t *mass2 = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass2, (vector_t){A/2, 0});
    scene_add_body(scene, mass2);

    create_spring(scene, K, mass1, mass2);

    // double initial_energy = spring_potential(mass1, mass2, K);
    for (int i = 0; i < STEPS; i++) {
        // Want to ensure the center is always 0
        double x1 = body_get_centroid(mass1).x;
        double x2 = body_get_centroid(mass2).x;
        assert(within(1e-8, x1 + x2, 0));
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

    DO_TEST(test_spring_center);
    DO_TEST(test_drag_slows);
    DO_TEST(test_fast_gravity_conservation);
    DO_TEST(test_spring_conservation);

    puts("student_forces_test PASS");
}
