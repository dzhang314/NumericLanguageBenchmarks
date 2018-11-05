#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

double squared_norm(const double *v, size_t n) {
    double result = 0.0;
    for (size_t i = 0; i < n; ++i) result += v[i] * v[i];
    return result;
}

double squared_distance(const double *v, const double *w, size_t n) {
    double result = 0.0;
    for (size_t i = 0; i < n; ++i) result += (v[i] - w[i]) * (v[i] - w[i]);
    return result;
}

void normalize_points(double *points, size_t num_points, size_t dim) {
    for (size_t i = 0; i < num_points; ++i) {
        double *const point = points + i * dim;
        const double inv_norm = 1.0 / sqrt(squared_norm(point, dim));
        for (size_t j = 0; j < dim; ++j) point[j] *= inv_norm;
    }
}

void normalize_forces(const double *points, double *forces,
                      size_t num_points, size_t dim) {
    for (size_t i = 0; i < num_points; ++i) {
        double dot = 0.0;
        for (size_t j = 0; j < dim; ++j)
            dot += points[i * dim + j] * forces[i * dim + j];
        for (size_t j = 0; j < dim; ++j)
            forces[i * dim + j] -= dot * points[i * dim + j];
    }
}

double coulomb_energy(const double *points, size_t num_points, size_t dim) {
    double energy = 0.0;
    for (size_t i = 1; i < num_points; ++i)
        for (size_t j = 0; j < i; ++j)
            energy += 1.0 / sqrt(squared_distance(
                points + i * dim, points + j * dim, dim));
    return energy;
}

double coulomb_forces(const double *points, double *forces,
                      size_t num_points, size_t dim) {
    double energy = 0.0;
    double displ[dim];
    memset(forces, 0, num_points * dim * sizeof(double));
    for (size_t i = 1; i < num_points; ++i) {
        for (size_t j = 0; j < i; ++j) {
            double dist_sq = 0.0;
            for (size_t k = 0; k < dim; ++k) {
                displ[k] = points[i * dim + k] - points[j * dim + k];
                dist_sq += displ[k] * displ[k];
            }
            const double inv_dist = 1.0 / sqrt(dist_sq);
            energy += inv_dist;
            const double inv_dist3 = inv_dist / dist_sq;
            for (size_t k = 0; k < dim; ++k) displ[k] *= inv_dist3;
            for (size_t k = 0; k < dim; ++k) forces[i * dim + k] += displ[k];
            for (size_t k = 0; k < dim; ++k) forces[j * dim + k] -= displ[k];
        }
    }
    return energy;
}

double quadratic_line_search(const double *points,
                             double energy,
                             const double *step_direction,
                             double step_size,
                             double *workspace,
                             size_t num_points, size_t dim) {
    for (size_t i = 0; i < num_points * dim; ++i)
        workspace[i] = points[i] + step_size * step_direction[i];
    normalize_points(workspace, num_points, dim);
    double new_energy = coulomb_energy(workspace, num_points, dim);
    if (new_energy < energy) {
        for (;;) {
            step_size *= 2.0;
            for (size_t i = 0; i < num_points * dim; ++i)
                workspace[i] = points[i] + step_size * step_direction[i];
            normalize_points(workspace, num_points, dim);
            const double newer_energy = coulomb_energy(
                workspace, num_points, dim);
            if (newer_energy >= new_energy) return 0.25 * step_size *
                (4.0 * new_energy - newer_energy - 3.0 * energy) /
                (2.0 * new_energy - newer_energy - energy);
            new_energy = newer_energy;
        }
    } else {
        for (;;) {
            step_size *= 0.5;
            for (size_t i = 0; i < num_points * dim; ++i)
                workspace[i] = points[i] + step_size * step_direction[i];
            normalize_points(workspace, num_points, dim);
            const double newer_energy = coulomb_energy(
                workspace, num_points, dim);
            if (newer_energy < energy) return 0.5 * step_size *
                (new_energy - 4.0 * newer_energy + 3.0 * energy) /
                (new_energy - 2.0 * newer_energy + energy);
            new_energy = newer_energy;
        }
    }
}

int main(int argc, char **argv) {

    if (argc < 3) {
        printf("Usage: %s <input-file-name> <num-steps>\n", argv[0]);
        return EXIT_FAILURE;
    }

    FILE *input_file = fopen(argv[1], "r");
    if (input_file == NULL) {
        printf("Error: could not open input file %s.\n", argv[1]);
        return EXIT_FAILURE;
    }

    size_t num_points, dim;
    fscanf(input_file, "%zu %zu", &num_points, &dim);
    double *const points = malloc(num_points * dim * sizeof(double));
    double *const forces = malloc(num_points * dim * sizeof(double));
    double *const workspace = malloc(2 * num_points * dim * sizeof(double));
    if (points == NULL || forces == NULL) {
        printf("Error: could not allocate enough memory.\n");
        return EXIT_FAILURE;
    }

    for (size_t i = 0; i < num_points * dim; ++i)
        fscanf(input_file, "%lf", points + i);

    const int num_steps = atoi(argv[2]);
    double step_size = 1.0e-3;
    double energy = coulomb_forces(points, forces, num_points, dim);
    normalize_forces(points, forces, num_points, dim);

    for (int step = 0; step < num_steps; ++step) {
        step_size = quadratic_line_search(points, energy, forces, step_size, 
                                          workspace, num_points, dim);
        for (size_t i = 0; i < num_points * dim; ++i)
            points[i] += step_size * forces[i];
        normalize_points(points, num_points, dim);
        energy = coulomb_forces(points, forces, num_points, dim);
        normalize_forces(points, forces, num_points, dim);
    }

    printf("%+.16e\n", energy);
    for (size_t i = 0; i < num_points; ++i) {
        for (size_t j = 0; j < dim; ++j) {
            printf("%+.16e", points[i * dim + j]);
            putchar(j + 1 == dim ? '\n' : '\t');
        }
    }

    return EXIT_SUCCESS;

}
