/* must be defined here */
#define NCLASS 3
#define NLGS   4

int nlocks_lg[] = {1,1,1,1} ;
int init_pos[] = {2, 0, 3};

/* only lock routing */
int routs[NCLASS][NLGS][2] = {
    {
	{0, 0},
	{0, 0},
	{3, 2},
	{2, 1},
    },
    {
	{1, 39},
	{2, 5},
	{0, 2},
	{0, 0},
    },
    {
	{0, 0},
	{0, 0},
	{0, 0},
	{3, 1}
    }
};

/* service times for locks AND local computation times */
int servs[NCLASS][2*NLGS] = {
    {
	0, 0, 0, 0, 30000, 7000000, 56000, 23000000
    },
    {
	7000,
	338000,
	6000,
	230000,
	30000,
	84000,
	0,
	0
    },
    {
	0,0,0,0,0,0,56000, 22000000
    }
};

double rout2[NCLASS][NLGS][NLGS] = {
    {{ 0.        ,  0.        ,  0.        ,  0.        },
     { 0.        ,  0.        ,  0.        ,  0.        },
     { 0.        ,  0.        ,  0.97701149,  0.02298851},
     { 0.        ,  0.        ,  0.01612903,  0.98387097}},
    {{ 0.98625833,  0.        ,  0.01374167,  0.        },
     { 0.        ,  0.88405797,  0.11594203,  0.        },
     { 0.24942966,  0.29201521,  0.45855513,  0.        },
     { 0.        ,  0.        ,  0.        ,  0.        }},
    {{ 0.,  0.,  0.,  0.},
     { 0.,  0.,  0.,  0.},
     { 0.,  0.,  0.,  0.},
     { 0.,  0.,  0.,  1.}}
};


/*
[[-- 338507.704986 --]
 [-- 7106.13352047 --]
 [-- 230605.906316 --]
 [-- 6378.42934783 --]
 [7204888.86874 83531.0878592 --]
 [28251.6777778 11492.3892668 --]
 [23009949.625 -- 22263582.6667]
 [56325.25 -- 34011.0]
 [4243584.58105 -- --]
 [12984.3650794 -- --]
 [-428.0 -- --]
 [15777.5 -- --]]


[array([[  0.,   0.,   0.,   0.,   0.,   0.],
       [  0.,   0.,   0.,   0.,   0.,   0.],
       [  0.,   0.,  85.,   3.,   2.,   0.],
       [  0.,   0.,   4.,   0.,   0.,   0.],
       [  0.,   0.,   1.,   0.,  61.,   0.],
       [  0.,   0.,   0.,   1.,   0.,   1.]]),
 array([[ 23541.,      0.,    328.,      0.,      0.,      0.],
       [     0.,   2928.,    384.,      0.,      0.,      0.],
       [   328.,    384.,    603.,      0.,      0.,      0.],
       [     0.,      0.,      0.,      0.,      0.,      0.],
       [     0.,      0.,      0.,      0.,      0.,      0.],
       [     0.,      0.,      0.,      0.,      0.,      0.]]),
 array([[ 0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  3.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.],
       [ 0.,  0.,  0.,  0.,  0.,  0.]])]                                
*/
