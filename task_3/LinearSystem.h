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

        // 0 = успех, отрицательные = ошибка
        int solve();
        void compute_residuals(double& r1, double& r2) const;

        bool naive_full_matrix_to_triangular();
        bool naive_traingular_solution();

        bool blocked_matrix_to_triangular();
        bool triangular_blocked_to_solution();
    private:
        bool init_matrix_from_formula(int s);
        bool init_matrix_from_file(char* file_name);
        bool init_perm();
        bool init_rhs();

        void get_block(int i, int j, double* dest) const;
        void set_block(int i, int j, const double* src);

        bool finding_block_pivot(int alpha, int& pivot_i, int& pivot_j, double* ws_block1, double* ws_block2, double* invert_block);

        void swap_blocked_rows(int i, int j);
        void swap_blocked_columns(int i, int j);

        double get_matrix_norm() const;
    private:
        double* A          = nullptr;
        double* b          = nullptr;
        double* solution   = nullptr;
        int*    perm       = nullptr; // перестановки матричных столбцов

        
        int* w_perm    = nullptr; // нужен для перестановок внутри обращаемого блока 
        double* block1 = nullptr;
        double* block2 = nullptr;
        double* block3 = nullptr;
        double* block4 = nullptr;

        int n = 0;
        int m = 0;
        int k = 0;
        int l = 0;
};

#endif
