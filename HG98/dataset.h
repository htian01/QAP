#pragma once
#include <stdlib.h>
#include <stdio.h>
#include <tuple>
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <cassert>
using namespace std;


tuple<float*, int> read_QAPLIB(const string file_name)
{
    ifstream file;
    int N = 0;
    string str;
    vector<string> v;
    string t;
    vector<int> W;
    vector<int> D;

    //read N
    file.open(file_name);
    getline(file, str);
    N = stoi(str);

    // read matrix
    while (getline(file, str))
    {
        istringstream in(str);
        while (in >> t) {
            v.push_back(t);
        }
    }
    int num_count = 0;
    for (string str : v)
    {
        if (num_count < N * N)
        {
            W.push_back(stoi(str));
        }
        else
        {
            D.push_back(stoi(str));
        }
        num_count++;
    }

    float* C = (float*)malloc(N * N * N * N * sizeof(float));
    for (int i = 0; i < N; i++)
    {
        for (int j = 0; j < N; j++)
        {
            for (int k = 0; k < N; k++)
            {
                for (int n = 0; n < N; n++)
                {
                    C[i * N * N * N + j * N * N + k * N + n] = float(W[i * N + k] * D[j * N + n]);
                }
            }
        }
    }

    return make_tuple(C, N);
}

class Model
{
public:
    Model(const int no_left, const int no_right, const int no_assignments, const int no_edges) : no_left(no_left), no_right(no_right), no_assignments(no_assignments), no_edges(no_edges) {}
    const int no_left;
    const int no_right;
    const int no_assignments;
    const int no_edges;
    vector<tuple<int, int, float>> assignments;
    vector<tuple<int, int, float>> edges;
    void add_assignment(const int id_assignment, const int id_left, const int id_right, const float cost)
    {
        assert(id_left < no_left);
        assert(id_right < no_right);
        assert(id_assignment == assignments.size());
        assignments.push_back(make_tuple(id_left, id_right, cost));
    }
    void add_edge(const int id_assignment1, const int id_assignment2, const float cost)
    {
        assert(edges.size() < no_edges);
        edges.push_back(make_tuple(id_assignment1, id_assignment2, cost));
    }
    void save_dd_format(const char* file_name, const float bound=0, const float primal=0)
    {
        assert(no_assignments == assignments.size());
        assert(no_edges == edges.size());
        FILE* fp = fopen(file_name, "w");
        fprintf(fp, "c %f %f\n", bound, primal);
        fprintf(fp, "p %d %d %d %d\n", no_left, no_right, no_assignments, no_edges);
        float cost;
        int id_left, id_right;

        for (int id_assignment = 0; id_assignment < assignments.size(); id_assignment++)
        {
            tie(id_left, id_right, cost) = assignments[id_assignment];
            fprintf(fp, "a %d %d %d %f\n", id_assignment, id_left, id_right, cost);
        }
        int id_assignment1, id_assignment2;

        for (auto it = edges.cbegin(); it != edges.cend(); it++)
        {
            tie(id_assignment1, id_assignment2, cost) = *it;
            fprintf(fp, "e %d %d %f\n", id_assignment1, id_assignment2, cost);
        }
        fclose(fp);
    }
};

Model* parse_dd_model(const char * filename)
{
    FILE* fp = fopen(filename, "r");
    if (!fp) { printf("Can't open %s\n", filename); exit(1); }
    int N[2], A, E, _A = 0, _E = 0;
    char LINE[256];
    float cost;
    Model* m = nullptr;
    while (fgets(LINE, sizeof(LINE) - 1, fp))
    {
        if (LINE[0] == 'p')
        {
            if (m || sscanf(&LINE[1], "%d %d %d %d\n", &N[0], &N[1], &A, &E) != 4) { printf("%s: wrong format1!\n", filename); exit(1); }
            m = new Model(N[0], N[1], A, E);
        }
        else if (LINE[0] == 'a')
        {
            int a, i0, i1;
            if (!m || sscanf(&LINE[1], "%d %d %d %f\n", &a, &i0, &i1, &cost) != 4
                || a != _A++ || _A > A || i0 < 0 || i0 >= N[0] || i1 < 0 || i1 >= N[1]) {
                printf("%s: wrong format2!\n", filename); exit(1);
            }
            m->add_assignment(a, i0, i1, cost);
        }
        else if (LINE[0] == 'e')
        {
            int a, b;
            if (!m || sscanf(&LINE[1], "%d %d %f\n", &a, &b, &cost) != 3
                || (_E++) >= E || a < 0 || a >= A || b < 0 || b >= A || a == b) {
                printf("%s: wrong format3!\n", filename); exit(1);
            }
            m->add_edge(a, b, cost);
        }
    }
    fclose(fp);
    return m;
}
Model* matrix_to_model(float* C, const int N)
{
    vector<tuple<int, int, float>> assignments;
    vector<tuple<int, int, float>> edges;
    float cost;
    for (size_t i = 0; i < N; i++)
    {
        for (size_t k = 0; k < N; k++)
        {
            assignments.push_back(make_tuple(i, k, C[i * N * N * N + k * N * N + i * N + k]));
        }
    }
    for (size_t i = 0; i < N; i++)
    {
        for (size_t j = i; j < N; j++)
        {
            for (size_t k = 0; k < N; k++)
            {
                for (size_t p = 0; p < N; p++)
                {
                    if (p != k)
                    {
                        cost = C[i * N * N * N + k * N * N + j * N + p] + C[j * N * N * N + p * N * N + i * N + k];
                        if (cost != 0)
                        {
                            edges.push_back(make_tuple(i * N + k, j * N + p, cost));
                        }
                    }
                }
            }
        }
    }
    Model*  m = new Model(N, N, assignments.size(), edges.size());
    m->assignments = assignments;
    m->edges = edges;
    return m;
}
tuple<float*, int> model_to_matrix(const Model* m)
{
    const int N = m->no_left;
    float* C = (float*)calloc(N * N * N * N, sizeof(float));
    int i, j, k, n;
    float cost;
    for (auto it = m->assignments.begin(); it != m->assignments.end(); it++)
    {
        tie(i, j, cost) = *it;
        C[i * N * N * N + j * N * N + i * N + j] = cost;
    }
    int a1, a2;
    int dummy;
    for (auto it = m->edges.begin(); it != m->edges.end(); it++)
    {
        tie(a1, a2, cost) = *it;
        tie(i, j, dummy) = m->assignments[a1];
        tie(k, n, dummy) = m->assignments[a2];
        C[i * N * N * N + j * N * N + k * N + n] = cost;
    }
    return make_tuple(C, N);
}
tuple<float*, int> read_dd_format(const string file_name)
{
    Model* m = parse_dd_model(file_name.c_str());
    return model_to_matrix(m);
}

