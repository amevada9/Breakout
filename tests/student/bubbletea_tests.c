#include <assert.h>
#include <math.h>
#include <stdlib.h>

#include "forces.h"
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

double spring_potential(double k, double x) {
    return 0.5 * k * x * x;
}

double kinetic_energy(body_t *body) {
    vector_t v = body_get_velocity(body);
    return body_get_mass(body) * vec_dot(v, v) / 2;
}

void test_energy_conservation_spring() {
    const double M = 10;
    const double K = 2;
    const double A = 3;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_centroid(mass, (vector_t){0, A});
    scene_add_body(scene, mass);
    body_t *anchor = body_init(make_shape(), INFINITY, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, anchor);
    create_spring(scene, K, mass, anchor);
    double initial_energy = spring_potential(K, A);
    for (int i = 0; i < STEPS; i++) {
        double position = A * cos(sqrt(K / M) * i * DT);
        double total_energy = spring_potential(K, position) + kinetic_energy(mass);
        // printf("Starting energy %f\n", initial_energy);
        // printf("Current energy %f\n", total_energy);
        assert(within(1e-4, total_energy / initial_energy, 1));
        scene_tick(scene, DT);
    }
    scene_free(scene);
}

// test the energy lost due to drag
void test_energy_lost_to_drag() {
    const double M = 5;
    const double V = 1000;
    const double GAMMA = 1.0;
    const double DT = 1e-6;
    const int STEPS = 1000000;
    scene_t *scene = scene_init();
    body_t *mass = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    body_set_velocity(mass, (vector_t){V, 0});
    scene_add_body(scene, mass);
    create_drag(scene, GAMMA, mass);
    double initial_energy = kinetic_energy(mass);
    double total_lost_energy = 0;
    // printf("Initial energy: %f\n", initial_energy);
    for (int i = 0; i < STEPS; i++) {
        // printf("Current kinetic energy: %f\n", kinetic_energy(mass));
        // printf("Total energy lost to drag: %f\n", total_lost_energy);
        // printf("Current energy + lost energy = %f\n", kinetic_energy(mass) +
        // total_lost_energy);
        assert(
            within(1e-4, (kinetic_energy(mass) + total_lost_energy) / initial_energy, 1));
        scene_tick(scene, DT);
        total_lost_energy +=
            GAMMA * body_get_velocity(mass).x * DT * body_get_velocity(mass).x;
    }
    scene_free(scene);
}

double gravity_potential(double G, body_t *body1, body_t *body2) {
    vector_t r = vec_subtract(body_get_centroid(body2), body_get_centroid(body1));
    return -G * body_get_mass(body1) * body_get_mass(body2) / sqrt(vec_dot(r, r));
}

// some other test that I need to write
void test_orbit() {
    const double m = 1, M = 1000000;
    const double G = 1e3;
    const double DT = 1e-6;
    const double R = 100;
    const int STEPS = 1000000;

    // populating scene with bodies and forces
    scene_t *scene = scene_init();
    body_t *satellite = body_init(make_shape(), m, (rgb_color_t){0, 0, 0});
    body_set_centroid(satellite, (vector_t){0, R});
    body_t *center = body_init(make_shape(), M, (rgb_color_t){0, 0, 0});
    scene_add_body(scene, center);
    vector_t v_orbit = (vector_t){sqrt(G * M / R), 0};
    body_set_velocity(satellite, v_orbit);
    scene_add_body(scene, satellite);
    create_newtonian_gravity(scene, G, satellite, center);
    double initial_energy = gravity_potential(G, satellite, center) + kinetic_energy(satellite);
    // printf("Initial energy = %f\n", initial_energy);
    bool was_period = false;
    double curr_period = 0;
    int period_ticks = 0;
    for (int i = 0; i < STEPS; i++) {
        /**
         * Center centroid test
         *
         * Check that the center body does not move much from its
         * original position at the origin.
         * The center will still move because the satellite causes a
         * gravitational force on the center
         */
        // printf("Center coords = (%f, %f)\n", body_get_centroid(center).x, body_get_centroid(center).y);
        // printf("Satellite coords = (%f, %f)\n", body_get_centroid(satellite).x, body_get_centroid(satellite).y);
        assert(vec_within(1e-2, body_get_centroid(center), (vector_t){0, 0}));

        /**
         * Satellite distance test
         *
         * Check that the distance between the center body and the satellite
         * stays close to its original distance (R = 100).
         * The distance will not be exactly constant because the center
         * moves slightly due to gravitaional force from the satellite.
        */
        vector_t r = vec_subtract(body_get_centroid(center), body_get_centroid(satellite));
        double T = 2 * M_PI * sqrt(pow(sqrt(vec_dot(r, r)), 3) / (G * M));
        // printf("Distance = %f\n", sqrt(vec_dot(r, r)));
        assert(within(1e-3, sqrt(vec_dot(r, r)) / R, 1));

        /**
         * Conservation of energy
         *
         * Energy should be conserved in this closed system of no outside forces.
         * Original energy consists of Gravitational energy and the kinetic
         * energy of the satellite.
         */
        double curr_energy = gravity_potential(G, satellite, center) + kinetic_energy(satellite) + kinetic_energy(center);
        // printf("Curr energy = %f\n", curr_energy);
        assert(within(1e-3, curr_energy / initial_energy, 1));

        /**
         * Orbital period test
         *
         * Period of orbit should obey Kepler's Third Law: T = 2 * pi * sqrt(r^3 / G * M)
         * Slight variation because of movement of center
         */
        if ((i * DT) > (T / 2) && body_get_centroid(satellite).y > 0 &&
            within(1e-2, body_get_centroid(satellite).x, body_get_centroid(center).x)) {
            was_period = true;
            period_ticks++;
        } else if (was_period) {
            was_period = false;
            // printf("Step: %d\n", i);
            // printf("Center coords = (%f, %f)\n", body_get_centroid(center).x,
            //        body_get_centroid(center).y);
            // printf("Satellite coords = (%f, %f)\n", body_get_centroid(satellite).x,
            //        body_get_centroid(satellite).y);
            // printf("Calculated period = %f\n", T);
            // printf("Actual period = %f\n", curr_period);
            // printf("\n");
            assert(within(1e-3, curr_period / T, 1));
            curr_period = (period_ticks / 2) * DT;
        }
        curr_period += DT;
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

    DO_TEST(test_energy_conservation_spring)
    DO_TEST(test_energy_lost_to_drag)
    DO_TEST(test_orbit)

    puts("student_test PASS");
}
