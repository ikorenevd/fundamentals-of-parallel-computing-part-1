#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H


// todo: добавить set для хранения статуса задачи
// нормальная обработка ошибок и более универсальных вывод
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

        // 0 = успех
        int solve();
        void compute_residuals(double& r1, double& r2) const;
    private:
        bool init_matrix_from_formula(int s);
        bool init_matrix_from_file(char* file_name);
        bool init_rhs();

        void get_block(int i, int j, double* dest) const;
        void set_block(int i, int j, const double* src);
    private:
        double* A = nullptr;
        double* b = nullptr;
        double* x = nullptr;
        double* norm_workspace = nullptr;

        int n = 0;
        int m = 0;
};

#endif
