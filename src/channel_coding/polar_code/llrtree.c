#include <stddef.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <stdbool.h>
#include <math.h>


#define CHECK_ARGS 0
#define INVALID_NAT UINT16_MAX
#define LOG_ACT 0
#define TRACE 0
#define COMPILE_AS_LIB 1
#if LOG_ACT
#define LOG(fmt, ...) do { \
    printf(fmt "\n", ##__VA_ARGS__); \
    fflush(stdout); \
} while(0)
#else
#define LOG(fmt, ...) ((void)0)
#endif
#if TRACE
    #include <execinfo.h>
    #include <signal.h>
#endif



typedef double real;
typedef u_int16_t nat;
typedef u_int8_t bit;
typedef u_int8_t (*dynamic_frozen_bit_at)(void* polarcode_ptr, bit* code_word, nat index, nat code_word_length);
#if TRACE
void handler(int sig) {
    void *array[50];
    size_t size;

    size = backtrace(array, 50);
    fprintf(stderr, "Signal %d:\n", sig);
    backtrace_symbols_fd(array, size, 2);
    exit(1);
}
#endif

#ifdef __cplusplus
extern "C" {
#endif

void* tree_list_init(nat n, nat paths_size);
void tree_list_deinit(void* vtreel);
void tree_list_reset(void* vtreel);
void tree_list_print(void* vlist);
void working_arg_print(void* vworking_arg);
int llrtreelist(
    real *init_llrs,
    void* vargs,
    void* vllr_tree,
    void* polarcode_ptr,
    bool *frozenset,
    bit *frozen_bits,
    bool *dynamic_frozen_set,
    dynamic_frozen_bit_at dynamic_frozen_bit_f
);
void* init_llrtreelist_working_args(nat n, nat L);
void deinit_llrtreelist_working_args(void* vargs);

#ifdef __cplusplus
}
#endif
/**
 * @brief returns log(exp(a) + exp(b))
 *
 * @param a a real number
 * @param b a real number
 *
 * @note this function is numerically stable a uses the formula
 * log(exp(a) + exp(b)) = max(a,b) + log(1 + exp(-|a-b|)))
 * For more on it, see, e.g.
 * https://en.wikipedia.org/wiki/LogSumExp
 *
 * @return log(exp(a) + exp(b))
 */
real log_sum_exp(real a, real b) {
    real max = fmax(a, b);
    real min = fmin(a, b);
    return max + log1p(exp(min - max));
}
/**
 * @brief return the min of a and b
 *
 * @param a a natural number
 * @param b a natural number
 *
 * @return the min of a and b
 */
static inline nat min(nat a, nat b) {
    return a >= b ? b : a;
}
/**
 * @brief returns the bit reverse of n over n_bits bits
 *
 * @param n the number to be bit reversed
 * @param n_bits the number of bit over which the reverse op is done
 *
 * @warning n should be such that 0 <= n < 2^{n_bits}. Otherwise, the 
 * output will be truncated. No argument checking is performed
 *
 * @return if n = b_{n_bits-1}···b_0, returns b_0···b_{n_bits-1}
 */
nat bit_rev(nat n, nat n_bits) {
    nat rev = 0;
    for (size_t i = 0; i < n_bits; i++) {
        rev <<= 1;
        rev |= n&1;
        n >>= 1;
    }
    return rev;
}
/**
 * @brief compute the metric associated with a given bit choice and a LLR
 *
 * @param u the chosen bit
 * @param llr the given llr
 *
 * @note For a given paths u = [u_0, u_1,..., u_i], the path metric defined as
 * PM(i) = metric(u_0, llr_0) + ··· + metric(u_i, llr_i), where 
 * llr_j = log(W_n^j(y,u[1:j-1]|0)) - log(W_n^j(y, u[1:j-1]|0)),
 * is such that for any two different paths u^l and u^k,
 * PM(i)^l > PM(i)^k if and only if W_n^i(y,u^l[1:i-1]|u_i^l) < W_n^i(y,u^k[1:i-1]|u_i^k)
 * Hence a smaller path metric means a better probability
 *
 * @return log(1 + exp(-(1-2u)llr))
 */
real metric(bit u, real llr) {
    real sign = u ? 1.0 : 0.0;
    return log_sum_exp(0, -(1-(2*sign))*llr);
}

/**
 * @brief a struct to hold an "index stack". It stores the available indexes, with pop and push 
 * operations. Since each index must be used at most once at the same time, this stack shouldn't
 * contain twice the same value. This is explicitely checked if CHECK_ARGS is set.
 */
typedef struct {
    size_t size;    /*The total size of the stack*/
    nat *stk;       /*The pointer to stack values*/
    ssize_t index;  /*The index at wich the stack points out*/
} index_stack_t;

/**
 * @brief just a type alias to handle the pointer more smoothly
 */
typedef index_stack_t* index_stack;
/**
 * @brief reset the stack to it's intial state, i.e filled
 * with [0, 1, 2, ..., stack->length-1]
 *
 * @param stack the stack to reinitialize
 *
 * @return 0 on success, -1 on error (if stack or stack->stk is NULL)
 */
int reset_index_stack(index_stack stack) {
#if CHECK_ARGS
    if (!stack || !stack->stk) return -1;
#endif
    //LOG("reset_index_stack : reseting index stack");
    for (size_t i = 0; i < stack->size; i++) {
        stack->stk[i] = i;
    }
    stack->index = 0;
    //LOG("   reset_index_stack : reseting index stack done");
    return 0;
}
/**
 * @brief init an index stack by allocating the required memory and data strucutre
 * 
 * @param stack_length the total length of the stack
 *
 * @return a index_stack which points to an index_stak_t structure. It is initialized as
 * [0, 1, 2, ..., stack_length-1]
 * The first poped value will be 0, then 1 and so on. Returns NULL if no memory is available
 *
 * @warning this function may return NULL. At the end of use of the stack, please use deinit_stack(). 
 */
index_stack init_index_stack(size_t stack_length) {
    nat *stk = calloc(stack_length, sizeof(nat));
    if (!stk) {
        LOG("Could not init stack");
        return NULL;
    }
    index_stack stack = calloc(1, sizeof(index_stack_t));
    if (!stack) {
        LOG("stack init : Could not init stack");
        free(stk);
        return NULL;
    }
    stack->index = 0;
    stack->size = stack_length;
    stack->stk = stk;
    reset_index_stack(stack);
    return stack;
}
/**
 * @brief pop the top value of the stack
 *
 * @param stack the stack to pop value
 *
 * @return the poped value on success. If the stack is empty, returns INVALID_NAT
 */
nat pop(index_stack stack) {
    ssize_t index = stack->index;
    size_t size = stack->size;
#if CHECK_ARGS
    if (index >= size) {
        LOG("stack pop : Stack is empty, cannot pop !!");
        return INVALID_NAT;
    }
#endif
    nat val = stack->stk[index];
    stack->stk[index] = INVALID_NAT;
    stack->index = index + 1;
    return val;
}
/**
 * @brief push a value on the top of the stack
 *
 * @param stck the stack to push value on
 * @param val the value to push on the stack
 *
 * @note If the stack is full, returns -1. If CHECK_ARGS is set, ensures
 * that the pushed value is not already present in the stack and returns
 * -1 if that is the case. 
 *
 * @return 0 on sucess, -1 if the stack is full and, if CHECK_ARGs is set, if 
 * the pushed value is already in the stack.
 */
int push(index_stack stack, nat val) {
    ssize_t index = stack->index;
    size_t size = stack->size;
#if CHECK_ARGS
    if (index < 1) {
        LOG("stack push : Stack is full, cannot push !!");
        return -1;
    }
    for (size_t i = stack->size - 1; i >= stack->index; i--) {
        if (stack->stk[i] == val) {
            LOG("stack push : cannot add %hu, it is already present");
            return -1;
        }
    }
#endif
    stack->stk[index - 1] = val;
    stack->index = index - 1;
    return 0;
}
/**
 * @brief deinit the allocated ressources of the stack
 *
 * @param stack the stack to deinit
 */
void deinit_stack(index_stack stack) {
    if (stack) {
        if (stack->stk) {
            free(stack->stk);
        }
        free(stack);
    }
    return;
}
/**
 * @brief a struct to hold the active path indexes. It is implemented as a simple 
 * fixed capacity array list. 
 */
typedef struct {
    size_t capacity;
    size_t size;
    nat *vals;
} active_index_l;
/**
 * @brief just a type alias
 */
typedef active_index_l* active_index;
/**
 * @brief reset the acitve index list to its initial state
 */
int reset_index_list(active_index list) {
#if CHECK_ARGS
    if (!list) return -1;
#endif
    //LOG("reset_index_list : reseting list size");
    list->size = 0;
    return 0;
}
/**
 * @brief init an empty active index list by allocating necessary ressources
 *
 * @param L the capacity of the list
 *
 * @return an active index list on success, NULL if no memory was available.
 * @warning after the list is no longer needed, please call deinit_active_index_list()
 */
active_index init_index_list(nat L) {
    nat *val = calloc(L, sizeof(nat));
    if (!val) {
        LOG("index list init : Could not init list");
        return NULL;
    }
    active_index index_list = calloc(1, sizeof(active_index_l));
    if (!index_list) {
        LOG("index list init : Could not init index list");
        free(val);
        return NULL;
    }
    index_list->capacity = L;
    index_list->size = 0;
    index_list->vals = val;
    return index_list;
}
/**
 * @brief deinit the active_index list by freeing the allocated ressources
 *
 * @param list the active_index list to free
 */
void deinit_active_index_list(active_index list) {
    if (!list) return;
    if (list->vals) {
        free(list->vals);
    }
    free(list);
    return;
}
/**
 * @brief add a new index to the active index list
 *
 * @param list the list to which add the active index to
 * @param val the value to add in the list
 *
 * @note if CHECK_ARGS is set, this function also checks
 * that val is not already present in the list an returns -1
 * if that's the case
 *
 * @return 0 on success, -1 on error, which can happen either
 * 1) if the index list is full
 * 2) if CHECK_ARGS is set and list is NULL or val is already in the list
 */
int add_active_index(active_index list, nat val) {
#if CHECK_ARGS
    if (!list) return -1;
    for (size_t i = 0; i < list->size; i++) {
        if (list->vals[i] == val) {
            LOG("index list add : Cannot add %hu, it is already in", val);
            return -1;
        }
    }
#endif
    if (list->size >= list->capacity) {
        LOG("index list add : List is full, cannot add new elements");
        return -1;
    }
    list->vals[list->size++] = val;
    return 0;
}
/**
 * @brief add all elements in the active index list
 *
 * @param list the list to add new indexes to
 * @param new_index a pointer to the indexes to add to list
 * @param new_index_size the number of new index to add to list
 *
 * @note if CHECK_ARGS is set, ensures that all elemetns can safely be
 * inserted and that no new elements is already in the list
 *
 * @return 0 on sucess, -1 on error, which can only happen if the arguments
 * are somewhat invalid and CHECK_ARGS is set
 */
int add_all_active_index(active_index list, nat *new_indexes, size_t new_indexes_size) {
#if CHECK_ARGS
    if (!list || !new_indexes) return -1;
    if (new_indexes_size > list->capacity) {
        LOG("List is full, cannot add new elements");
        return -1;
    }
    for (size_t i = 0; i < list->size; i++) {
        for (size_t j = 0; j < new_indexes_size; j++) {
            if (list->vals[i] == new_indexes[j]) {
                LOG("index list add all : cannot add %hu, it is already in", new_indexes[j]);
                return -1;
            }
        }
    }
#endif
    nat size = list->size;
    nat new_size = size + new_indexes_size;
    nat *vals = list->vals;

    for (size_t i = size; i < new_size; i++) {
        vals[i] = new_indexes[i - size];
    }
    list->size = new_size;
    return 0;
}
/**
 * @brief remove val from list
 *
 * @param list the list to remove val from
 * @param val the value to remove from list
 *
 * @note this function shift all elements in the array list
 * which are after val to the left. If iterating through the
 * list while removing elements, the iterations must be done 
 * in descending order, or otherwise elements won't be considered
 *
 * @return 0 on success, -1 if the element was not in the list
 */
int remove_active_index(active_index list, nat val) {
#if CHECK_ARGS
    if (!list) return -1;
#endif
    nat size = list->size;
    for (size_t i = 0; i < size; i++){
        if (list->vals[i] == val) {
           for (size_t j = i+1; j < size; j++) {
            list->vals[j-1] = list->vals[j];
           }
           list->size--;
           return 0;
        }
    }
    LOG("Couldn't remove index because it is not in list");
    return -1;
}
/**
 * @brief init a bit pointer of pointer of dim (n1, n2)
 *
 * @return a bit pointer of pointer on success, NULL if not
 * all the n1·n2 bits could not be allocated. 
 *
 * @note This function shouldn't be used to create a matrix as a 
 * mathematical object, as it creates fragmentation. It is designed
 * to consider each pointer as a separated array. 
 */
bit** init_bit_mat(nat n1, nat n2) {
    bit **mat = calloc(n1, sizeof(bit*));
    if (!mat) return NULL;
    for (size_t i = 0; i < n1; i++) {
        bit *row = calloc(n2, sizeof(bit));
        if (!row) {
            for (size_t j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
        mat[i] = row;
    }
    return mat;
}
/**
 * @brief init a bool pointer of pointer of dim (n1, n2)
 *
 * @return a bool pointer of pointer on success, NULL if not
 * all the n1·n2 bool could not be allocated. 
 *
 * @note This function shouldn't be used to create a matrix as a 
 * mathematical object, as it creates fragmentation. It is designed
 * to consider each pointer as a separated array. 
 */
bool** init_bool_mat(nat n1, nat n2) {
    bool **mat = calloc(n1, sizeof(bool*));
    if (!mat) return NULL;
    for (size_t i = 0; i < n1; i++) {
        bool *row = calloc(n2, sizeof(bool));
        if (!row) {
            for (size_t j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
        mat[i] = row;
    }
    return mat;
}
/**
 * @brief init a real pointer of pointer of dim (n1, n2)
 *
 * @return a real pointer of pointer on success, NULL if not
 * all the n1·n2 real could not be allocated. 
 *
 * @note This function shouldn't be used to create a matrix as a 
 * mathematical object, as it creates fragmentation. It is designed
 * to consider each pointer as a separated array. 
 */
real** init_real_mat(nat n1, nat n2) {
    real **mat = calloc(n1, sizeof(real*));
    if (!mat) return NULL;
    for (size_t i = 0; i < n1; i++) {
        real *row = calloc(n2, sizeof(real));
        if (!row) {
            for (size_t j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
        mat[i] = row;
    }
    return mat;
}
/**
 * @brief deinit a pointer of pointer of dim (n1, ?)
 */
void deinit_mat(void** mat, nat n1) {
    for (size_t i = 0; i < n1; i++) {
        if (mat[i]) free(mat[i]);
    }
    free(mat);
}
/**
 * @brief init an array or array of array pointers which are
 * inited to NULL. 
 *
 * @return an array or array (n1, n2) of array pointers which are
 * inited to NULL on success, NULL if not all the memory could 
 * be allocated.
 *
 * @warning it is the caller responsability to further init the
 * array pointers, but the function deinit_realarray_mat will 
 * free them.
 */
real*** init_realarray_mat(nat n1, nat n2) {
    real ***mat = calloc(n1, sizeof(real**));
    if (!mat) return NULL;
    for (size_t i = 0; i < n1; i++) {
        real **row = calloc(n2, sizeof(real*));
        if (!row) {
            for (size_t j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
        mat[i] = row;
    }
    return mat;
}

bool*** init_boolarray_mat(nat n1, nat n2) {
    bool ***mat = calloc(n1, sizeof(bool**));
    if (!mat) return NULL;
    for (size_t i = 0; i < n1; i++) {
        bool **row = calloc(n2, sizeof(bool*));
        if (!row) {
            for (size_t j = 0; j < i; j++) free(mat[j]);
            free(mat);
            return NULL;
        }
        mat[i] = row;
    }
    return mat;
}
void deinit_boolarray_mat(bool*** mat, nat n1, nat n2) {
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2; j++) {
            //LOG("deinit_realarray_mat : Freeing mat[i][j] = %zu", mat[i][j]);
            if (mat[i][j]) free(mat[i][j]);
        }
        //LOG("deinit_realarray_mat : Freeing mat[i] = %zu", mat[i]);
        free(mat[i]);
    }
    //LOG("deinit_realarray_mat : Freeing mat = %zu", mat);
    free(mat);
}
/**
 * @brief deinit an array or array of array pointers of dim (n1, n2)
 *
 * @warning the arrays at mat[n1-1][n2-1] are not freeed. It is the
 * caller responsability to deallocate them. However, every other allocated
 * array is freed. 
 */
void deinit_realarray_mat(real*** mat, nat n1, nat n2) {
    for (size_t i = 0; i < n1; i++) {
        for (size_t j = 0; j < n2-1; j++) {
            //LOG("deinit_realarray_mat : Freeing mat[i][j] = %zu", mat[i][j]);
            if (mat[i][j]) free(mat[i][j]);
        }
        //LOG("deinit_realarray_mat : Freeing mat[i] = %zu", mat[i]);
        free(mat[i]);
    }
    //LOG("deinit_realarray_mat : Freeing mat = %zu", mat);
    free(mat);
}
/**
 * @brief a struct to hold the llr tree list final results
 */
typedef struct {
    bit **paths;    /*(path_size, n) The final paths*/
    real **llrs;    /*(path_size, n) The final llrs associated with each path*/
    real *metrics;  /*(path_size) The final paths metrics associated with each path*/
    nat n;          /*The size of each path*/
    nat paths_size; /*The number of selected paths*/
} llr_tree_list;

typedef llr_tree_list * tree_list;

void tree_list_print(void* vlist) {
    if (!vlist) return;
    tree_list list = (tree_list) vlist;
    printf("llrs :");
    for (nat l = 0; l < list->paths_size; l++) {
        printf("    path : %d with metric %f | ", l, list->metrics[l]);
        for (nat i = 0; i < list->n; i++) {
            printf("%f ", list->llrs[l][i]);
        }
        printf("\n");
    }
}

void* tree_list_init(nat n, nat paths_size) {
    bit **final_paths = init_bit_mat(paths_size, n);
    if (!final_paths) return NULL;
    real *final_metrics = calloc(paths_size, sizeof(real));
    if (!final_metrics) {
        deinit_mat((void**)final_paths,paths_size);
        return NULL;
    }
    real **final_llrs = init_real_mat(paths_size, n);
    if (!final_llrs) {
        deinit_mat((void**)final_paths,paths_size);
        free(final_metrics);
        return NULL;
    }
    tree_list res = calloc(1, sizeof(llr_tree_list));
    if (!res) {
        deinit_mat((void**)final_llrs,paths_size);
        deinit_mat((void**)final_paths,paths_size);
        free(final_metrics);
        return NULL;
    }
    res->llrs = final_llrs;
    res->metrics = final_metrics;
    res->paths = final_paths;
    res->paths_size = paths_size;
    res->n = n;

    return res;
}
void tree_list_reset(void* vtreel) {
    tree_list treel = (tree_list) vtreel;
    if (!treel) return;
    bit** paths = treel->paths;
    real** llrs = treel->llrs;
    if (paths != NULL && llrs != NULL) {
        for(nat l = 0; l < treel->paths_size; l++) {
            memset(paths[l], 0, treel->n*sizeof(bit));
            memset(llrs[l], 0, treel->n*sizeof(real));
        }
    }
    real* metrics = treel->metrics;
    if (metrics) memset(metrics, 0, treel->paths_size*sizeof(real));
    return;
}
void tree_list_deinit(void* vtreel) {
    tree_list treel = (tree_list) vtreel;
    if (!treel) return;
    if (treel->llrs) deinit_mat((void**)treel->llrs,treel->paths_size);
    if (treel->paths) deinit_mat((void**)treel->paths,treel->paths_size);
    if (treel->metrics) free(treel->metrics);
    free(treel);
    return;
}
/**
 * @brief deinit a llr tree list
 */
void deinit_llrlist(llr_tree_list treelist) {
    if (treelist.llrs) {
        deinit_mat((void**) treelist.llrs, treelist.paths_size);
    }
    if (treelist.paths) {
        deinit_mat((void**) treelist.paths, treelist.paths_size);
    }
    if (treelist.metrics) {
        free(treelist.metrics);
    }
    return;
}
/**
 * The following two functions could be stacked (but what's the risk of stack overflow ?)
 */
/**
 * @brief return the even indexed bit of a vecotr indexed by [1:u_length]
 *
 * @param u the bit vector to extact the even-indexed bit with first index 1
 * @param u_length the length of u
 *
 * @return a bit pointer containing the even indexed bit of u, NULL on memory 
 * failure
 *
 * @warning the even index are considered in the indexing [1:u_length]. In c, this
 * is translated to [0:u_length-1]. Hence, the even indexes are in c the odd ones.
 */
bit* even(bit* u, nat u_length) {
#if CHECK_ARGS
    if (!u) return NULL;
#endif
    nat new_length = u_length / 2;
    bit* new_u = calloc(new_length, sizeof(bit));
    if (!new_u) return NULL;

    size_t j = 0;
    //LOG("even : doing computation for length %hu", u_length);
    for (size_t i = 1; i < u_length; i+=2) {
        new_u[j++] = u[i];
    }
    return new_u;
}
bit* even_buff(bit* u, nat u_length, bit* dst) {
#if CHECK_ARGS
    if (!u || !dst) return NULL;
#endif
    nat new_length = u_length / 2;
    size_t j = 0;
    //LOG("even : doing computation for length %hu", u_length);
    for (size_t i = 1; i < u_length; i+=2) {
        *dst++ = u[i];
    }
    return dst;
}
/**
 * @brief given u, returns even(u) ^ odd(u)
 *
 * @param u the bit vecotr to xor
 * @param u_lenght the length of u
 *
 * @return even(u) ^ odd(u) = [u_0 ^ u_1, u_2 ^ u_3, ...]
 *
 * @note if u_length is odd, the last bit is ignored
 */
bit* xor(bit* u, nat u_length) {
#if CHECK_ARGS
    if (!u) return NULL;
#endif
    nat new_length = u_length / 2;
    bit* new_u = calloc(new_length, sizeof(bit));
    if (!new_u) return NULL;
    //LOG("xor : doing cumputations for length %hu", u_length);
    size_t j = 0;
    for (size_t i = 1; i < u_length; i+=2) {
        new_u[j++] = u[i-1] ^ u[i];
    }

    return new_u;
}
bit* xor_buffer(bit* u, nat u_length, bit* dst) {
#if CHECK_ARGS
    if (!u || !dst) return NULL;
#endif
    nat new_length = u_length / 2;
    for (size_t i = 1; i < u_length; i+=2) {
        *dst++ = u[i-1] ^ u[i];
    }

    return dst;
}
/**
 * @brief compute llr_mat[l][d][i] using the polar recursion.
 *
 * @param bit the bit vecotr on which the llrs are dependent
 * @param len_bits the length of bit
 * @param i the index of which to compute the llr, @warning in bit reversed !
 * @param d the depth of the recursion. Please set it to 0 for a complete calculation
 * @param l the path index of bit
 * @param llr_mat the llr table of arrays overwhich computation is done
 * @param n the dim of the arrays llr_mat contains pointer of
 * @param L the list length, it is such that 0 <= l < L.
 *
 * @returns 0 on success, -1 on failure
 */
int compute(
    bit *bits,
    size_t len_bits,
    nat i,
    nat d,
    nat l,
    real ***llr_mat,
    bool ***llr_in,
    nat n,
    nat L,
    bit *subvec_buffer
) {
#if CHECK_ARGS
    if (!bits || !llr_mat) return -1;
    if (!((0 <= i) && (i <= n))) {
        LOG("compute : Invalid i=%hu against n=%hu", i, n);
        return -1;
    }
#endif
#if LOG_ACT
    if (d != 0) {
        //LOG("compute : ============ Compute for i=%hu, d=%hu, n=%hu ==========", i, d, n);
    } else {
       // LOG("BEGIN-RECURSION =================== i=%hu, d=%hu, n=%hu ==========", i, d, n);
    }
#endif
    nat shift = n / (1 << (d+1));
    //LOG("compute : Computed shift(n=%hu,d=%hu) = %hu", n,d, shift);
    nat new_length = len_bits / 2;
    int ret = 0;
    bit *subvec = NULL;
    nat raw_i = d != 0 ? i / (n >> d) : i;
    //LOG("compute : Computed raw_i(i=%hu, n=%hu, d=%hu) = %hu", i, n, d, raw_i);
    //LOG("compute : Getting array and path %hu and depth %hu for i=%hu", l, d, i);
    real *arr = llr_mat[l][d];

    bool *b_arr = llr_in[l][d];
    //LOG("compute : Cheking if present at %hu", i);
    if (!b_arr[raw_i]) { 
        //stack the sub u ?
        subvec = xor_buffer(bits, len_bits, subvec_buffer);
        if (!subvec) {
            //LOG("compute : Couldn't allocate space for a xored vector");
            return -1;
        }
        ret = compute(
            subvec_buffer,
            new_length,
            i,
            d+1,
            l,
            llr_mat,
            llr_in,
            n,
            L,
            subvec);
        //LOG("compute : ============= end of call =================");
        subvec = NULL;
        //LOG("compute : From i=%hu, computed new i=%hu at depth=%hu", i, new_i, d);
        subvec = even_buff(bits, len_bits, subvec_buffer);
        if (!subvec) {
            //LOG("compute : Couldn't allocate space for a even vector");
            return -1;
        }
        nat new_i = i + ((i % (2*shift)) < shift ? shift : -shift);
        //LOG("compute : Computed new_i(i=%hu,shift=%hu) = %hu", i, shift, new_i);
        ret = compute(
            subvec_buffer,
            new_length,
            new_i,
            d+1,
            l,
            llr_mat, 
            llr_in,
            n,
            L,
            subvec);
        //LOG("compute : ============= end of call =================");
        subvec = NULL;
        if (ret < 0) return -1;
        //LOG("compute : Getting computed array and path %hu and depth %hu", l, d+1);
        real *other_array = llr_mat[l][d+1];
        nat index = i / (n >> (d+1));
        //LOG("compute : Computed index(i=%hu, n=%hu, d=%hu) = %hu", i, n, d, index);
        real l2 = other_array[index];
        nat new_index = new_i / (n >> (d+1));
        //LOG("compute : Computed new_index(new_i=%hu, n=%hu, d=%hu) = %hu", new_i, n, d, new_index);
        real l1 = other_array[new_index];
        if (len_bits % 2 == 0) {
            arr[raw_i] = log_sum_exp(l1+l2, 0) - log_sum_exp(l1, l2);
        } else {
            arr[raw_i] = bits[len_bits-1] == 0 ? (l2 + l1) : (l2 - l1);
        }

        b_arr[raw_i] = true;
        if (index != n-1 && index == ((1 << (d+1)) -1)) {
            //LOG("compute : For i=%hu, d=%hu, index=%hu  : reset upper array", i, d, index);
            memset(llr_in[l][d+1], 0, (1 << (d+1)));
        } else {
           // LOG("compute : For i=%hu, d=%hu, index=%hu  : do NOT reset upper array", i, d, index);
        }
    }
    //LOG("comptue : retrun %f", llr_mat[l][d][raw_i]);
    return 0;
}
/**
 * @brief duplicate path at index l, update metrics and llr_mat
 *
 * @param paths the array of arrays pointer of each path
 * @param l the path index to duplicate
 * @param i the bit index to write to path[l] and path[new_path_index]
 * @param metrics the metric table keeping track of the goodness of each active path
 * @param n the maximal length of each path
 * @param L the list size, it is such that 0 <= l < L
 * @param stack the stack to pop a new available path index
 * @param new_metrcis an array of size 2 which contains the new metrics for the paths
 *
 * @return the new path index on success, -1 on failure
 */
int duplicate_paths(
    bit **paths,
    nat l,
    nat i,
    real *metrics,
    nat n,
    real ***llr_mat,
    bool ***llr_in,
    nat depth,
    nat L,
    index_stack stack,
    real *new_metrics 
    ) {
#if CHECK_ARGS
    if (!paths || !metrics || !llr_mat || !stack) return -1;
#endif
    nat new_path_index = pop(stack);
    //LOG("duplicate_paths : retrieved new index=%hu", new_path_index);
    if (new_path_index == INVALID_NAT) return -1;
    //LOG("duplicate_paths : fetching and copying new path");
    bit *new_path = paths[new_path_index];
    bit *old_path = paths[l];
    memcpy(new_path, old_path, i*sizeof(bit));
    //LOG("duplicate_paths : updating boths paths bits and metrics");
    new_path[i] = 0;
    old_path[i] = 1;
    metrics[new_path_index] = new_metrics[0];
    metrics[l] = new_metrics[1];

    //LOG("duplicate_paths : copying path precise values (%hu->%hu)", l, new_path_index);
    //TODO : if iterated on bitrev order, only need to copy i entries!!
    memcpy(llr_mat[new_path_index][0], llr_mat[l][0], n*sizeof(real));
    memcpy(llr_in[new_path_index][0], llr_in[l][0], n*sizeof(bool));
    for (nat d = 1; d < depth-1; d++) {
        real *arr = llr_mat[new_path_index][d];
        memcpy(arr, llr_mat[l][d], (1 << d)*sizeof(real));

        bool *b_arr = llr_in[new_path_index][d];
        memcpy(b_arr, llr_in[l][d], (1 << d)*sizeof(bool));
    }
    return new_path_index;
}
/**
 * @brief a comparison function to sort a table or reals
 */
int comp(const void *a, const void *b) {
    real da = *(const real*)a;
    real db = *(const real*)b;
    return (da > db) - (da < db);
}
/**
 * @brief set active_metrics[i] = true for the best thresholds new_metrics
 *
 * @param new_metrics the new computed metrics, contains invalid values
 * @param active_metrcis indicate whether new_metrics[i] is a valid one
 * @param size length of new_metrics and active_metrics
 * @param threshold the number of new metrics to retain
 *
 * @note this function is in averge O(size·log(size)) (depends on the used libc and qsort). 
 * It could be done in O(size) by using the median of median method. But for our practical purposes, 
 * size is relatively small (in [1, 256]) and hence this doesn't make such a difference. 
 */
int select_treshold_best_metrics(real *new_metrics, bool* active_metrics, nat size, nat threshold) {
#if CHECK_ARGS
    if (!new_metrics || !active_metrics) return -1;
    if (!((0 <= threshold) && (threshold <= size))) {
        LOG("select_treshold_best_metrics : invalid threshold, cannot remove all paths");
        return -1;
    }
#endif
    real copy[size];
    nat copy_index = 0;
    for (nat i = 0; i < size; i++){
        if (active_metrics[i]) {
            copy[copy_index++] = new_metrics[i];
        }
    }
    real t;
    nat selected = 0;
    if (copy_index >= threshold) {
        qsort(copy, copy_index, sizeof(real), &comp);
        t = copy[threshold-1];
    } else {
        //LOG("select_treshold_best_metrics : for threshold %hu and active metrics number %hu, no work is needed", threshold, copy_index);
        return 0;
    }
    //LOG("select_treshold_best_metrics : computed numerical threshold is %f", t);
    for (nat i = 0; i < size; i++) {
        if (active_metrics[i]) {
            if (selected < threshold && new_metrics[i] <= t) {
                //LOG("select_treshold_best_metrics : secect index %hu with metric %f", i, new_metrics[i]);
                selected++;
            } else {
                //LOG("select_treshold_best_metrics : remove index %hu with metric %f", i, new_metrics[i]);
                active_metrics[i] = false;
            }
        }
    }
    //LOG("select_treshold_best_metrics : out of %hu, selected %hu paths (threshold was %hu)", size, selected, threshold);
    return 0;
}
/**
 * @brief a struct to hold llrtreelist algorithm working buffer and reduce the number of allocations of them
*/
typedef struct {
    real ***llr_mat;
    bool ***llr_in;
    bit **paths;
    real *metrics;
    void* stack;
    void* index_list;
    real *new_metrics;
    bool *active_new_metrics;
    bit* subvec_buffer;
    nat L;
    nat n;
    nat depth;
} llr_tree_list_working_struct_t;
typedef llr_tree_list_working_struct_t* working_arg;

void working_arg_print(void* vworking_arg) {
    if (!vworking_arg) return;
    working_arg args = (working_arg) vworking_arg;
    for (nat l = 0; l < args->L; l++) {
        for (nat d = 0; d < args->depth; d++) {
            if (!args->llr_mat[l][d]) continue;
            for (nat i = 0; i < args->n; i++) {
                printf("%f ", args->llr_mat[l][d][i]);
            }
            printf("\n");
        }
    }
}

bool working_arg_is_inited(working_arg args) {
    if (!args) return false;
    return !(!args->llr_mat || !args->llr_in || !args->paths 
        || !args->metrics || !args->stack || !args->index_list || !args->new_metrics
        || !args->active_new_metrics || !args->subvec_buffer);
}

void deinit_llrtreelist_working_args(void* vargs) {
    if (!vargs) return;
    working_arg args = (working_arg) vargs;
    if (args->llr_mat) {
        deinit_realarray_mat(args->llr_mat, args->L, args->depth);
        args->llr_mat = NULL;
    }
    if (args->llr_in) {
        deinit_boolarray_mat(args->llr_in, args->L, args->depth);
        args->llr_in = NULL;
    }
    if (args->paths) {
        deinit_mat((void**) args->paths, args->L);
        args->paths = NULL;
    }
    if (args->metrics) {
        free(args->metrics);
        args->metrics = NULL;
    }
    if (args->stack) {
        deinit_stack(args->stack);
        args->stack = NULL;
    }
    if (args->index_list) {
        deinit_active_index_list(args->index_list);
        args->index_list = NULL;
    }
    if (args->new_metrics) {
        free(args->new_metrics);
        args->new_metrics = NULL;
    }
    if (args->active_new_metrics) {
        free(args->active_new_metrics);
        args->active_new_metrics = NULL;
    }
    if (args->subvec_buffer) {
        free(args->subvec_buffer);
        args->subvec_buffer = NULL;
    }
    args->L = 0;
    args->n = 0;
    args->depth = 0;
    free(args);
    return;
}
void* init_llrtreelist_working_args(nat n, nat L){

    working_arg res = calloc(1, sizeof(llr_tree_list_working_struct_t));

    nat depth = (nat) log2((real)n) + 1;
    res->llr_mat = init_realarray_mat(L, depth);
    if (!res->llr_mat) goto err;
    res->llr_in = init_boolarray_mat(L, depth);
    if (!res->llr_in) goto err;
    res->paths = init_bit_mat(L, n);
    if (!res->paths) goto err;
    res->metrics = calloc(L, sizeof(real));
    if (!res->metrics) goto err;
    res->stack = init_index_stack(L);
    if (!res->stack) goto err;
    res->index_list = init_index_list(L);
    if (!res->index_list) goto err;
    res->new_metrics = calloc(2*L, sizeof(real));
    if (!res->new_metrics) goto err;
    res->active_new_metrics = calloc(2*L, sizeof(bool));
    if (!res->active_new_metrics) goto err;
    res->subvec_buffer = calloc(n, sizeof(bit));
    if (!res->subvec_buffer) goto err;

    res->L = L;
    res->depth = depth;
    res->n = n;

    for (nat l = 0; l < L; l++) {
        res->llr_mat[l][0] = calloc(n, sizeof(real));
        if (!res->llr_mat[l][0]) goto err;

        res->llr_in[l][0] = calloc(n, sizeof(bool));
        if (!res->llr_in[l][0]) goto err;
        for (nat d = 1; d < depth-1; d++){
            res->llr_mat[l][d] = calloc((1 << d), sizeof(real));
            if (!res->llr_mat[l][d]) goto err;

            res->llr_in[l][d] = calloc((1 << d), sizeof(bool));
            if (!res->llr_in[l][d]) goto err;
        }

        res->llr_in[l][depth-1] = calloc(n, sizeof(bool));
        if (!res->llr_in[l][depth-1]) goto err;
        memset(res->llr_in[l][depth-1], 1, n*sizeof(bool));
    }

    return res;
err:
    deinit_llrtreelist_working_args(res);
    return NULL;
}
int reset_working_args(working_arg args, real *init_llrs) {
#if CHECK_ARGS
    if (!init_llrs || !working_arg_is_inited(args)) return -1;
#endif
    //LOG("reset_working_args : on reseting ");
    for (nat l = 0; l < args->L; l++) {
        //LOG("reset_working_args : reseting llr_in at length %hu", l);
        memset(args->llr_in[l][0], 0, args->n*sizeof(bool));
        for (nat d = 1; d < args->depth-1; d++){
           // LOG("   reset_working_args : reseting llr_in at depth %hu", d);
            memset(args->llr_in[l][d], 0, (1 << d)*sizeof(bool));
        }
        //the last array should never be modified...
        //memset(args->llr_in[l][args->depth-1], 1, args->n*sizeof(bool));
        args->llr_mat[l][args->depth-1] = init_llrs;
    }
    memset(args->metrics, 0, args->L*sizeof(real));
    reset_index_list(args->index_list);
    reset_index_stack(args->stack);
    return 0;
}
/**
 * @brief compute an llrtree, that is the L best paths, with their associated metrics and llrs
 *
 * @param init_lrrs the initial values of the likelihood ratios, i.e log(P(X=0)/P(X=0)) + log(W(y_i|0)/W(y_i|1))
 * @param n the blocklength, i.e len(init_llrs)
 * @param L the list size. The SCL decoding will greedily consider the L most likely paths
 * @param frozenset the set of frozen indices. i in a frozen index iff frozenset[i] = true
 * @param frozenbit the value of the frozen bits (padded with zeros to reduce overhead)
 *
 * @return a vlid llrtree on success, a zero inited one on failure   
*/
int llrtreelist(
    real *init_llrs,
    void* vargs,
    void* vllr_tree,
    void* polarcode_ptr,
    bool *frozenset,
    bit *frozen_bits,
    bool *dynamic_frozen_set,
    dynamic_frozen_bit_at dynamic_frozen_bit_f
) {
    working_arg args = (working_arg) vargs;
    tree_list llr_tree = (tree_list) vllr_tree;
#if CHECK_ARGS
    if (!init_llrs || !frozenset || !frozen_bits || !working_arg_is_inited(args) || !vllr_tree) goto err;
    if (args->n != (1 << (int)log2((real)args->n))) goto err;
    if (!polarcode_ptr || !dynamic_frozen_bit_f || !dynamic_frozen_set) goto err;
#endif


    reset_working_args(args, init_llrs);
    tree_list_reset(vllr_tree);

    nat depth = args->depth;
    nat n = args->n;
    nat L = args->L;

    real*** llr_mat = args->llr_mat;
    bool*** llr_in = args->llr_in;
    bit** paths = args->paths;
    real* metrics = args->metrics;
    index_stack stack = args->stack;
    active_index index_list = args->index_list;
    real* new_metrics = args->new_metrics;
    bool* active_new_metrics = args->active_new_metrics;
    bit *subvec_buffer = args->subvec_buffer;

#if LOG_ACT
    LOG("INPUT:");
    printf("llrs|");
    for (nat i = 0; i < n; i++) {
        printf("%f ",init_llrs[i]);
    }
    printf("\n");
    printf("metrics|");
    for (nat i = 0; i < L; i++) {
        printf("%f ", args->metrics[i]);
    }
    printf("\n");
#endif
    //LOG("llrtreelist : starting to work");
    nat first_path_i = pop(stack);
    if (first_path_i == INVALID_NAT) goto err;
    if (add_active_index(index_list,first_path_i) < 0) goto err;

    int ret = 0;
    for (nat i = 0; i < n; i++) {
       // LOG("llrtreelist : Into main loop for i=%d", i);
        nat rev = bit_rev(i, depth - 1);
        //LOG("   llrtreelist : Computed rev index : %hu", rev);
        if (frozenset[i] || dynamic_frozen_set[i]) {
            //LOG("llrtreelist : Froren bit detected");
            bit ui = frozen_bits[i];
            LOG("   llrtreelist : Writting frozen bit %d at depth %hu", ui, i);
            for (nat l_index = 0; l_index < index_list->size; l_index++) {
                nat l = index_list->vals[l_index];
                //LOG("llrtreelist : Computing LLR for l=%d, i=%d", l, i);
                ret = compute(
                    paths[l], 
                    i, 
                    rev,
                    0, 
                    l, 
                    llr_mat,
                    llr_in,
                    n,
                    L,
                    subvec_buffer);
                if (ret < 0) goto err;
                LOG("       llrtreelist : accesing path %hu", l);
                ui = dynamic_frozen_set[i] ? dynamic_frozen_bit_f(polarcode_ptr, paths[l], i+1, n) : ui;
                paths[l][i] = ui;
                LOG("       llrtreelist : accesing metrics %hu", l);
                metrics[l] += metric(ui, llr_mat[l][0][rev]);
            }
            LOG("   llrtreelist : Writting frozen bit %d at depth %hu done", ui, i);
        } 
        else {
            LOG("llrtreelist : Compute the new metrics");
            //compute the new metrics...
            nat selectn = min(2*index_list->size, L);
            for (nat l_index = 0; l_index < index_list->size; l_index++) {
                nat l = index_list->vals[l_index];
                //LOG("llrtreelist : Computing LLR for l=%d, i=%d", l, i);
                    ret = compute(
                    paths[l], 
                    i, 
                    rev, 
                    0,
                    l, 
                    llr_mat,
                    llr_in,
                    n,
                    L,
                    subvec_buffer);
                if (ret < 0) goto err;
                real llr = llr_mat[l][0][rev];
                real old_metric = metrics[l];
                for (bit ui = 0; ui < 2; ui++) {
                    LOG("llrtreelist : Duplicating metrics");
                    real new_met = old_metric + metric(ui, llr); 
                    new_metrics[2*l + ui] = new_met;
                    active_new_metrics[2*l + ui] = true;
                }
            }
            //duplicate all paths;
            if (2*index_list->size <= L) {
                LOG("llrtreelist : Duplicate all paths...From l=%zu", index_list->size);
                nat new_paths_indexes[index_list->size];
                for (nat l_index = 0; l_index < index_list->size; l_index++) {
                    nat l = index_list->vals[l_index];
                    ret = duplicate_paths(
                        paths,
                        l,
                        i,
                        metrics,
                        n,
                        llr_mat,
                        llr_in,
                        depth,
                        L,
                        stack,
                        new_metrics + 2*l
                    );
                    if (ret < 0) goto err;
                    new_paths_indexes[l_index] = ret;
                }
                ret = add_all_active_index(index_list, new_paths_indexes, index_list->size);
                if (ret < 0) goto err;
                LOG("llrtreelist : New used index list size is %zu", index_list->size);
            } else {
                LOG("llrtreelist : Selecting good paths. From l=%zu", index_list->size);
                if (select_treshold_best_metrics(new_metrics, active_new_metrics, 2*L, selectn) < 0) goto err;
                for (ssize_t l_index = (ssize_t)index_list->size - 1; l_index >= 0; l_index--) {
                    nat l = index_list->vals[l_index];
                    if (!active_new_metrics[2*l] && !active_new_metrics[2*l + 1]) {
                        //LOG("   llrtreelist : deleting path %hu", l);
                        if (remove_active_index(index_list, l) < 0) goto err;
                        if (push(stack, l) < 0) goto err;
                    }
                }
                for (ssize_t l_index = (ssize_t)index_list->size - 1; l_index >= 0; l_index--) {
                    nat l = index_list->vals[l_index];
                    if (active_new_metrics[2*l] && !active_new_metrics[2*l + 1]) {
                       // LOG("   llrtreelist : updating on branch path %hu, choose bit 0", l);
                        paths[l][i] = 0;
                        metrics[l] = new_metrics[2*l];
                    } else if (!active_new_metrics[2*l] && active_new_metrics[2*l + 1]) {
                        //LOG("   llrtreelist : updating on branch path %hu, choose bit 1", l);
                        paths[l][i] = 1;
                        metrics[l] = new_metrics[2*l + 1];
                    } else {
                        LOG("   llrtreelist : duplicating path %hu", l);
                        ret = duplicate_paths(
                            paths,
                            l,
                            i,
                            metrics,
                            n,
                            llr_mat,
                            llr_in,
                            depth,
                            L,
                            stack,
                            new_metrics + 2*l
                        );
                        if (ret < 0) goto err;
                        if (add_active_index(index_list, ret) < 0) goto err;
                    }
                }
                LOG("llrtreelist : to l=%zu", index_list->size);
            }
            memset(active_new_metrics, 0, 2*L);
        }
    }
    LOG("llrtreelist : init final values");
    real **final_llrs = llr_tree->llrs;
    bit **final_paths = llr_tree->paths;
    real *final_metrics = llr_tree->metrics;

    //FIXME : in the future, it might happen that not all paths are allocated 
    LOG("llrtreelist : filling final values");
    for (nat l_index = 0; l_index < index_list->size; l_index++) {
       // LOG("llrtreelist : getting path index for path-index-index %hu", l_index);
        nat l = index_list->vals[l_index];
       // LOG("llrtreelist : copying metrics %hu", l);
        final_metrics[l_index] = metrics[l];
        real *llrs = final_llrs[l_index];
        //LOG("llrtreelist : copying llrs and path : %hu", l);
        for (nat i = 0; i < n; i++){
            //LOG("llrtreelist : copying llrs for i=%hu and irev=%hu", i, bit_rev(i, depth-1));
            final_paths[l_index][i] = paths[l][i];
            llrs[i] = llr_mat[l][0][bit_rev(i, depth - 1)];
        }
    }
    LOG("llrtreelist : filling unused inices");
    for (nat l_index = index_list->size; l_index < L; l_index++) {
        final_metrics[l_index] = 0;
    }
   // LOG("llrtreelist : sucess");
#if LOG_ACT
    LOG("DECODED :");
    for (nat l = 0; l < args->L; l++) {
        printf("path:%hu|", l);
        for (nat i = 0; i < n; i++) {
            printf("%hu", llr_tree->paths[l][i]);
        }
        printf("\nllr:%hu|", l);
        for (nat i = 0; i < n; i++) {
            printf("%f ", llr_tree->llrs[l][i]);
        }
        printf("\nmetric:%f\n", llr_tree->metrics[l]);
    }
#endif
    return 0;
err:
    return -1;
}
#if !COMPILE_AS_LIB
int main(void){
#if TRACE
    signal(SIGSEGV, handler);
    signal(SIGABRT, handler);
#endif
    LOG("Here we starts uwwu");
    nat n = 128;
    nat l = 16;
    real init_llrs[n];
    bool frozenset[n];
    bit frozen_bits[n];
    for (size_t i = 0; i < n; i++) {
        long rdm = random();
        frozenset[i] = (rdm >> 2)&1;
        frozen_bits[i] = rdm&1;
        real val = ((real)(rdm & 32) / (real)32) + 1;
        init_llrs[i] = val; 
    }
    LOG("Here we to compute !!");
    llr_tree_list llrtree = llrtreelist(init_llrs, n, l,frozenset,frozen_bits);
    LOG("We're done");
    if (llrtree.llrs) {
        deinit_mat((void**) llrtree.llrs, llrtree.paths_size);
    }
    if (llrtree.paths) {
        deinit_mat((void**) llrtree.paths, llrtree.paths_size);
    }
    if (llrtree.metrics) {
        free(llrtree.metrics);
    }
    return 0;
}
#endif