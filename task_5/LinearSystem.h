#ifndef LINEAR_SYSTEM_H
#define LINEAR_SYSTEM_H


// todo: добавить set для хранения статуса задачи
// нормальная обработка ошибок и более универсальных вывод
class LinearSystem
{
    public:
        LinearSystem(int _matrix_size, int _block_size);
        ~LinearSystem();

        // 1 - если успешно, 0 - если нет, ничего не инициализируем и памяти не касаемся
        bool memory_alloc();
        void free_memory();

        bool init(int s, char* file_name);

        void print_matrix(int r) const;
        void print_rhs(int r) const;
        void print_solution(int r) const;

        // 0 = успех, отрицательные = ошибка
        int solve();
        void compute_residuals(double& r1, double& r2) const;

        // 1 - если успешно, 0 - если нет
        bool blocked_matrix_to_triangular();
        bool triangular_blocked_to_solution();
    private:
        // 1 - если успешно, 0 - если нет
        bool init_matrix_from_formula(int s);
        // 1 - если успешно, 0 - если нет
        bool init_matrix_from_file(char* file_name);
        // 1 - если успешно, 0 - если нет
        bool init_perm();
        // 1 - если успешно, 0 - если нет
        bool init_rhs();

        void get_block(int i, int j, double* dest) const;
        void set_block(int i, int j, const double* src);

        bool finding_block_pivot(int alpha, int& pivot_i, int& pivot_j, double* ws_block1, double* ws_block2, double* invert_block);

        void swap_blocked_rows(int i, int j, int first_column);
        void swap_blocked_columns(int i, int j);

        double get_matrix_norm() const;

        double* get_block_address(int i, int j) const;
    private:
        double* coeff_matrix    = 0;
        double* right_hand_side = 0;
        double* solution        = 0;
        int*    perm            = 0; // перестановки матричных столбцов
        
        int* w_perm    = nullptr; // нужен для перестановок внутри обращаемого блока 
        double* block1 = nullptr;
        double* block2 = nullptr;
        double* block3 = nullptr;

        int matrix_size      = 0;
        int block_size       = 0;
        int full_block_count = 0;
        int remainder_size   = 0;
        int block_count      = 0;

        double original_matrix_norm = 0.;
};

#endif
