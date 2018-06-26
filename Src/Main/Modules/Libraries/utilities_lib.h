void test_start(char *file,int line);
int countf(void *p);
void init_func_PtoV(sFunc_PtoV_T ***func);
void add_func_PtoV(sFunc_PtoV_T ***func,void (*f)(Patch_T *patch),char *task,Coord_T coord);
void run_func_PtoV(sFunc_PtoV_T **func,char *task,Patch_T *patch);
Coord_T find_coord(Patch_T *patch);
Collocation_T get_collocation(char *coll);

#define TEST_START test_start(__FILE__,__LINE__);


void IJK(int l, int *n, int *i, int *j, int *k);
int L(int *n, int i, int j, int k);
int I(int l, int *n);
int J(int l, int *n);
int K(int l, int *n);
int IsThisEdge(int *n,int p);
