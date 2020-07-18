def Kronecker_Product(A, B): # A otimes B
    m1, n1 = len(A), len(A[0])
    m2, n2 = len(B), len(B[0])
    # res       = [[A[0][0]*B, A[0][1]*B, ..., A[0][n1-1]*B], 
    #              [A[1][0]*B, A[1][1]*B, ..., A[1][n1-1]*B],
    #              ...,
    #              [A[m1-1][0]*B, A[m1-1][1]*B, ..., A[m1-1][n1-1]*B]]
    #
    # res[i][j] = A[i//m2][j//n2] * B[i%m2][j%n2]
    res = [[0 for _ in range(n1*n2)] for __ in range(m1*m2)]
    for i in range(m1*m2):
        for j in range(n1*n2):
            res[i][j] = A[i//m2][j//n2] * B[i%m2][j%n2]
    return res

def linear_equation_solver(A, b): # find vector x s.t. Ax = b
    m, n = len(A), len(A[0])
    ext_mat = [[0 for _ in range(n+1)] for __ in range(m)]
    for i in range(m):
        for j in range(n):
            ext_mat[i][j] = A[i][j]
        ext_mat[i][n] = b[i]
    pos = 0
    for i in range(n+1):
        nonzero = pos
        flg = True
        while flg:
            if nonzero >= m:
                flg = False
                break
            else:
                if ext_mat[nonzero][i] != 0:
                    break
                else:
                    nonzero += 1
        if flg:
            keep = ext_mat[nonzero][i]
            for j in range(i, n+1):
                ext_mat[nonzero][j], ext_mat[pos][j] = ext_mat[pos][j], ext_mat[nonzero][j]
                ext_mat[pos][j] /= keep
            ext_mat[pos][i] = 1
            for j in range(pos+1, m):
                raw_keep = ext_mat[j][i]
                for k in range(i, n+1):
                    ext_mat[j][k] -= raw_keep * ext_mat[pos][k]
            pos += 1
    exist_solution = True
    for i in range(m):
        sub_flg = False
        if ext_mat[i][-1]:
            for j in range(n):
                if ext_mat[i][j] != 0:
                    sub_flg = True
        else:
            sub_flg = True
        exist_solution *= sub_flg
    if exist_solution:
        for j in range(m, 0, -1):
            flg = True
            pos = 0
            while flg:
                if pos > n:
                    flg = False
                    break
                else:
                    if ext_mat[j-1][pos] != 0:
                        break
                    else:
                        pos += 1
            if flg:
                for k in range(j-1):
                    raw_keep = ext_mat[k][pos]
                    for l in range(pos, n+1):
                        ext_mat[k][l] -= raw_keep * ext_mat[j-1][l]
        solution = {
            "solution": [0 for _ in range(n)], 
            "subspace": []
        }
        pointer = [-1 for _ in range(n)]
        cnt = 0
        for i in range(m):
            for j in range(n):
                if ext_mat[i][j] == 1:
                    cnt += 1
                    solution["solution"][j] = ext_mat[i][n]
                    pointer[j] = i
                    break
        solution["subspace"] = [[0 for _ in range(n)] for _ in range(n-cnt)]
        D = {}
        pnt = 0
        for i in range(n):
            if pointer[i] == -1:
                solution["subspace"][pnt][i] = 1
                D[pnt] = i
                pnt += 1
        for i in range(n):
            if pointer[i] != -1:
                for j in range(n-cnt):
                    solution["subspace"][j][i] = -ext_mat[pointer[i]][D[j]]
        return solution
    else:
        return {"solution": [], "sub_solution": []}

"""N, M = map(int, input().split())
b = list(map(int, input().split())) # len(b) = N
A = [list(map(int, input().split())) for _ in range(N)]
print(linear_equation_solver(A, b))"""
