import numpy as np
import re
import matplotlib.pyplot as plt


def dot2pctable(dotfile):
    parent_child_table = []
    file = open(dotfile, encoding='utf8')
    line = file.readline().strip()
    if len(line) > 1:
        line = file.readline().strip()
        match = re.search(r'\d+\s->\s\d+.*', line)
        if match != None:
            pc_pair = re.match(r'(\d+)\s->\s(\d+).*(\d)\.',line).groups()
            parent_child_table.append(np.array(list(pc_pair)).astype(int))
    return np.array(parent_child_table)

def one_hot_encoding(index_list, category_list):
    encoding = np.zeros((len(category_list), len(index_list)))
    encoding[:, index_list] = 1
    return encoding

def snv_assign(C_CNV, Q, A, E, U, F):
    """
    the function for assigning unsampled SNVs to the trees, using brutal force with minimum
    distance criteria to identify the possible branch and allele of a SNV given
    n - number of clones
    m - number of samples
    l - number of SVs
    g - number of sampled SNVs
    g_un - number of unsampled SNVs
    r - number of CNVs
    :param C_CNV: n*2r allelic specific CNV
    :param Q: g_un * r mapping matrix which maps the unsampled SNVs to CNV segments, q_ij=1 if ith SNV maps to jth CNV
    :param A: n*n, a_ij = 1 if i is the ancestor of j, diagonal is 0, which means i is not the ancestor of i
    :param U: m*n frequency matrix
    :param F: m*g_un frequency matrix
    :return:
    """
    n, r = C_CNV.shape
    g_un = Q.shape[0]
    r = int(r/2)
    clone_idx_range = range(0, n)
    C_hat_1 = np.dot(C_CNV[:, :r], np.transpose(Q)) # n*g_un, the copy number of CNV at SNV position
    C_hat_2 = np.dot(C_CNV[:, r:], np.transpose(Q))  # n*g_un, the copy number of CNV at SNV position
    C_hat_1_parent = np.dot(E.T, C_hat_1)
    C_hat_2_parent = np.dot(E.T, C_hat_2)
    min_dist = np.full((g_un), np.inf)
    min_node = np.full((g_un), -1)
    for b in clone_idx_range:
        print(b)
        ### normal copy number=1
        C_SNV_clone_1 = C_hat_1[b, :] # g_un
        C_SNV_clone_2 = C_hat_2[b, :]
        C_SNV_clone_parent_1 = C_hat_1_parent[b, :]
        C_SNV_clone_parent_2 = C_hat_2_parent[b, :]
        valid_snv_idx = np.array(list(set(np.append(np.where(C_SNV_clone_1 == 1)[0],np.where(C_SNV_clone_1 - C_SNV_clone_parent_1 > 1)[0]))))
        F_est = U[:,b] * C_SNV_clone_1[valid_snv_idx] + np.dot(U, A[b, :][:,np.newaxis]* C_hat_1[:,valid_snv_idx])
        #print('1',F_est, F[:, valid_snv_idx])
        dist = np.sum(np.abs(F_est - F[:, valid_snv_idx]),axis=0)
        dist_stack = np.column_stack((min_dist[valid_snv_idx], dist))
        argmin = np.argmin(dist_stack, axis=-1)
        if (argmin == 1).any():
            print(min_node[valid_snv_idx[argmin == 1]])
            min_node[valid_snv_idx[argmin == 1]] = b
            print(min_node[valid_snv_idx[argmin == 1]])
            min_dist[valid_snv_idx] = np.min(dist_stack, axis=-1)
            print(min_node, min_dist)

        valid_snv_idx = np.array(list(set(np.append(np.where(C_SNV_clone_2 == 1)[0],np.where(C_SNV_clone_2 - C_SNV_clone_parent_2 > 1)[0]))))
        F_est = U[:, b] * C_SNV_clone_2[valid_snv_idx] + np.dot(U, A[b, :][:, np.newaxis] * C_hat_2[:, valid_snv_idx])
        #print('2',F_est, F[:, valid_snv_idx])
        dist = np.sum(np.abs(F_est - F[:, valid_snv_idx]),axis=0)
        dist_stack = np.column_stack((min_dist[valid_snv_idx], dist))
        argmin = np.argmin(dist_stack, axis=-1)
        if (argmin == 1).any():
            print(min_node[valid_snv_idx[argmin == 1]])
            min_node[valid_snv_idx[argmin == 1]] = b
            print(min_node[valid_snv_idx[argmin == 1]])
            min_dist[valid_snv_idx] = np.min(dist_stack, axis=-1)
            print(min_node, min_dist)

        ### copy number > 1
        valid_snv_idx = np.where(C_SNV_clone_1 > 1)[0]
        F_est = U[:, b][:,np.newaxis] + np.dot(U, A[b, :][:, np.newaxis] * C_hat_1[:, valid_snv_idx] / C_SNV_clone_1[valid_snv_idx])
        #print('3',F_est, F[:, valid_snv_idx])
        dist = np.sum(np.abs(F_est - F[:, valid_snv_idx]),axis=0)
        dist_stack = np.column_stack((min_dist[valid_snv_idx], dist))
        argmin = np.argmin(dist_stack, axis=-1)
        if (argmin == 1).any():
            print(min_node[valid_snv_idx[argmin == 1]])
            min_node[valid_snv_idx[argmin == 1]] = b
            print(min_node[valid_snv_idx[argmin == 1]])
            min_dist[valid_snv_idx] = np.min(dist_stack, axis=-1)
            print(min_node, min_dist)

        valid_snv_idx = np.where(C_SNV_clone_2 > 1)[0]
        F_est = U[:, b][:,np.newaxis] + np.dot(U, A[b, :][:, np.newaxis] * C_hat_2[:, valid_snv_idx] / C_SNV_clone_2[valid_snv_idx])
        #print('4',F_est, F[:, valid_snv_idx])
        dist = np.sum(np.abs(F_est - F[:, valid_snv_idx]),axis=0)
        dist_stack = np.column_stack((min_dist[valid_snv_idx], dist))
        argmin = np.argmin(dist_stack, axis=-1)
        if (argmin == 1).any():
            print(min_node[valid_snv_idx[argmin == 1]])
            min_node[valid_snv_idx[argmin == 1]] = b
            print(min_node[valid_snv_idx[argmin == 1]])
            min_dist[valid_snv_idx] = np.min(dist_stack, axis=-1)
            print(min_node, min_dist)

    return min_node, min_dist



if __name__ == '__main__':
    ### test case
    n=3
    m=1
    g_un=5
    k=4
    np.random.seed(0)
    U = np.array([[0.1, 0.5, 0.3, 0.1]])
    C_CNV = np.array([[1,2,1,1,1,1,4,1],
                      [1,2,1,3,2,2,1,1],
                      [1,2,1,1,1,1,1,1],
                      [1,1,1,1,1,1,1,1],])
    A = np.array([[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [1, 1, 0, 0],
                  [1, 1, 1, 0]])
    Q = np.eye(4)
    C_SNV = np.array([[1, 1, 4, 0],
                      [0, 1, 1, 1],
                      [0, 1, 1, 0],
                      [0, 0, 0, 0]])
    E = np.array([[0, 0, 0, 0],
                  [0, 0, 0, 0],
                  [1, 1, 0, 0],
                  [0, 0, 1, 0]])
    F = np.dot(U, C_SNV)
    min_node, min_dist = snv_assign(C_CNV, Q, A, E, U, F)
