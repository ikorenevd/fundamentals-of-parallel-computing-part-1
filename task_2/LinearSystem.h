#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H

class LinearSystem
{
    public:
        LinearSystem(int n, int m);
        ~LinearSystem();

        bool memory_alloc();
        void free_memory();

        bool init(int s, char* file_name);

        void print_matrix(int r) const;
        void print_rhs(int r) const;
        void print_solution(int r) const;

        bool solve();
        void compute_residuals(double& r1, double& r2) const;
    private:
        bool init_matrix_from_formula(int s);
        bool init_matrix_from_file(char* file_name);
        bool init_rhs();
        bool init_solution();

        void get_block(int i, int j, double* dest) const; 

        void calculate_matrix_norm();
    private:
        double* A = nullptr;
        double* b = nullptr;
        double* x = nullptr;

        int n = 0;
        int m = 0;

        double matrix_norm = 0.;    
};

#endif