from gurobipy import *
import math
import numpy as np
import time
import networkx as nx
import random

def createGraph(input_file):
    global G, n, m, a0, start_time, apsp_time, n_nodes, n_edges, conn_comp, av_degree, density, card_comp
    
    G = nx.Graph()
    for j in range(0,n):
        G.add_node(j)
   
    f = open(input_file, "r")
    string = f.readline()
    for i in range(0, m):
        string = f.readline()
        string = string.split()
        j = int(string[0])-1
        k = int(string[1])-1
        G.add_edge(j, k)
    f.close()
    
    for i in range(0, n):
        G.add_edge(i, i)
    
    print("graph created")
                
    conn_comp = nx.number_connected_components(G)
    n_nodes   = len(nx.nodes(G))
    n_edges   = len(nx.edges(G))
    av_degree = 0
    for e in nx.degree(G):
        av_degree = av_degree + e[1]
    av_degree = av_degree / n_nodes
    density   = nx.density(G)
    print("conn. comp.: " + str(conn_comp))
    print("num vertices: " + str(n_nodes))
    print("num edges: " + str(n_edges))
    print("density: " + str(density))
    print("av degree: " + str(av_degree))
    print("---Compute all pairs shortest path - running time: %s seconds ---" % (time.time() - start_time))

def run():
    global total_runtime, runtime, n, m, a0, feasible, best_sequence, I_input, n_ff, T_input
    try:
        m = Model("mip1")
        m.Params.outputFlag = 1  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 2); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        
        # ------------------------------INPUT --------------------------------
        I = I_input  # Fire starts from vertex 0
        b0 = []
        d0 = []
        for i in range(n):
            if i in I:
                b0.append(1) # These are the vertices burned at time j=0 ----------------------------(10)
            else:
                b0.append(0)
            d0.append(0) # No vertex is defended at time j=0 ---------------------------------------(11)
        a0 = [] # Graph's adjacency matrix ---------------------------------------------------------(12)
        for i in range(n):
            list = []
            for j in range(n):
                list.append(0)
            a0.append(list)
        for i in range(n):
            for j in range(n):
                if i in G.neighbors(j):
                    a0[i][j] = 1
        
        T = T_input
        
        #---------------------------- VARIABLES -----------------------------------------------------------
        
        b = []   #-----------------------------------------------------------------------------------------(A:b)
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            b.append(temp)
        for i in range(n):
            for j in range(T):
                b[i][j] = m.addVar(vtype=GRB.BINARY, name="b,%s" % str(i+1) + "," + str(j+1))
        
        d = []   #-----------------------------------------------------------------------------------------(A:d)
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            d.append(temp)
        for i in range(n):
            for j in range(T):
                d[i][j] = m.addVar(vtype=GRB.BINARY, name="d,%s" % str(i+1) + "," + str(j+1))
        
        #a = []   #-----------------------------------------------------------------------------------------(B)
        #for j in range(n):
        #    temp1 = []
        #    for i in range(n):
        #        temp2 = []
        #        for k in range(n):
        #            temp2.append(0)
        #        temp1.append(temp2)
        #    a.append(temp1)
        #for j in range(n):
        #    for i in range(n):
        #        for k in range(n):
        #            a[j][i][k] = m.addVar(vtype=GRB.BINARY, name="a,%s" % str(j+1) + "," + str(i+1) + "," + str(k+1))
        
        #---------------------------- CONSTRAINTS ---------------------------------------------------------
        
        for i in range(n): #-------------------------------------------------------------------------------(2)
            for j in range(T):
                if j == 0:
                    m.addConstr(b[i][j] >= b0[i])
                else:
                    m.addConstr(b[i][j] >= b[i][j-1])
        
        for i in range(n): #-------------------------------------------------------------------------------(3)
            for j in range(T):
                if j == 0:
                    m.addConstr(d[i][j] >= d0[i])
                else:
                    m.addConstr(d[i][j] >= d[i][j-1])
                    
        for i in range(n): #-------------------------------------------------------------------------------(4 y 5)
            for j in range(T):
                for k in range(n):
                    if j == 0:
                        m.addConstr(b[i][j] + d[i][j] >= b0[k] * a0[i][k]) 
                    else:
                        m.addConstr(b[i][j] + d[i][j] >= b[k][j-1] * a0[i][k])
                        m.addConstr(b[i][j] + d[i][j] >= b[k][j-1] * a0[i][k])
                        #m.addConstr(b[i][j] + d[i][j] >= b[k][j-1] * a0[i][k] * (1-d[i][j-1]))
                        #m.addConstr(b[i][j] + d[i][j] >= b[k][j-1] * a0[i][k] * (1-d[k][j-1]))                        
                         
        for i in range(n): #-------------------------------------------------------------------------------(7)
            for j in range(T):
                if j == 0:
                    m.addConstr(b0[i] + d0[i] <= 1)
                else:
                    m.addConstr(b[i][j] + d[i][j] <= 1)
                    
        for i in range(n): #-------------------------------------------------------------------------------(7 y 8)
            for j in range(T):
                sum_ = 0
                if j == 0:
                    for k in range(n):
                        sum_ = sum_ + b0[k] * a0[i][k]
                    m.addConstr(sum_ >= b[i][j])
                else:
                    for k in range(n):
                        sum_ = sum_ + b[k][j-1] * a0[i][k]
                        #sum_ = sum_ + b[k][j-1] * a0[i][k] * (1-d[i][j-1])
                        #sum_ = sum_ + b[k][j-1] * a0[i][k] * (1-d[k][j-1])
                    m.addConstr(sum_ >= b[i][j])
                    
        for j in range(T): #-------------------------------------------------------------------------------(9)
            if j == 0:
                sum_ = 0
                for i in range(n):
                    sum_ = sum_ + d[i][j] - d0[i]
                m.addConstr(sum_ <= nff)                
            else:
                sum_ = 0
                for i in range(n):
                    sum_ = sum_ + d[i][j] - d[i][j-1]
                m.addConstr(sum_ <= nff)
                
        # --------------------------
        # Se va a necesitar una matriz d para cada firefighter
        #for i in range(n):
        #    for j in range(1,T):
        #        sum_ = 0
        #        for k in range(n):
        #            if j == 1:
        #                sum_ = sum_ + a0[i][k] * d[k][j-1]
        #            else:
        #                sum_ = sum_ + a0[i][k] * (d[k][j-1] - d[k][j-2])
        #                #sum_ = sum_ + a0[i][k] * (d[k][j-1])
        #        m.addConstr(d[i][j] <= sum_ + d[i][j-1])
                        
        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------
        
        b_transpose = np.array(b).T.tolist()
        m.setObjective(sum(b_transpose[T-1]), GRB.MINIMIZE)#-----------------------------------------------(1)
                
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        #m.computeIIS()
        #m.write("model_or.lp")
        
        m.optimize()
        runtime = m.Runtime
        print("The run time is %f" % runtime)
        if m.status == GRB.INFEASIBLE:
            #print('Model is infeasible')
            feasible = False
        else:
            print("Obj:", m.objVal)
            feasible = True
            b_out = []
            d_out = []
            sol_d = []
            sol_b = []
            a_out = []
            burned_vertices = 0
            for v in m.getVars():
                varName = v.varName
                varNameSplit = varName.split(',')
                if varName[0] == 'd' and varNameSplit[2] == str(n):
                   sol_d.append(v.x)
                if varName[0] == 'b' and varNameSplit[2] == str(n):
                   sol_b.append(v.x)
                if varName[0] == 'b':
                    b_out.append(v.x)
                if varName[0] == 'd':
                    d_out.append(v.x)
                if varName[0] == 'a':
                    a_out.append(v.x)
            
            # Get burned vertices and
            # sequence of firefighters
            d = []
            for i in range(n):
                temp = []
                for j in range(T):
                    temp.append(0)
                d.append(temp)                
            j = 0
            k = 0
            #><
            for i in range(len(d_out)):
                if d_out[i] > 0.9:
                    d[j][k] = 1
                else:
                    d[j][k] = 0
                if (i+1)%T == 0:
                    j = j + 1
                    k = -1
                k =  k + 1
                

            # Get burned vertices and
            # sequence of firefighters
            b = []
            for i in range(n):
                temp = []
                for j in range(T):
                    temp.append(0)
                b.append(temp)
            j = 0
            k = 0
            for i in range(len(b_out)):
                if b_out[i] > 0.9:
                    b[j][k] = 1
                else:
                    b[j][k] = 0
                if (i+1)%T == 0:
                    j = j + 1
                    k = -1
                k =  k + 1
                             
            #for i in range(n):
            #    print()
            #    for j in range(n):
            #        print(a[i][j])
            
            b_transpose = np.array(b).T.tolist()
            not_change = True
            i = 0
            while(not_change):
                if i == 0:
                    if sum(b0) == sum(b_transpose[i]):
                        not_change = False
                    else: 
                        i += 1     
                else:
                    if sum(b_transpose[i]) == sum(b_transpose[i-1]):
                        not_change = False
                    else:
                        i += 1
            if i == 0: # Este es un caso atípico en el que solo se quema un nodo
                i = 1
                
            # Si ya no hay más nodos quemados en la columna (tiempo) j es por una de dos:
            # 1) Ya no hay más nodos que quemar.
            # 2) Existe una siquiente ronda en la que se defienden más nodos
            #    en el vecindario del fuego.
            # Lo siguiente detecta el caso 2)
            flag2 = 0
            for j in range(n):
                # Checar si algún nodo recién defendido
                if d[j][i] == 1 and d[j][i-1] == 0:
                    # es vecino de alguien con fuego
                    for p in range(n):
                        if a0[j][p] * (1-d[j][i-1]) == 1 and a0[j][p] * (1-d[j][i-1]) == 1  \
                                                        and b[p][i-1] == 1:
                       # if a[i-1][j][p] == 1 and b[p][i-1] == 1:
                            flag2 = 1
            if flag2 == 1:
                i = i + 1
                
            print()
            n_rounds = i
            print("Number of rounds: " + str(n_rounds))
            print("Burned vertices: " + str(sum(b_transpose[T-1])))
            
            d_transpose = np.array(d).T.tolist()
            ff_seq = []
            for i in range(n_rounds):
                temp = []
                for j in range(n):
                    if d_transpose[i][j] == 1:
                        temp.append(j)
                ff_seq.append(temp)
            #print(ff_seq)
            
            out_seq = []
            out_seq.append(set(ff_seq[0]))
            for i in range(1, len(ff_seq)):
                s1 = set(ff_seq[i-1])
                s2 = set(ff_seq[i])
                s3 = s2 - s1
                out_seq.append(s3)
            print("Firefighters sequence: " + str(out_seq))
            
            # --------------- DEBUGGING ------------------------
            
            debug = False # To debug, set to True
            
            if debug == True:
                print()
                print("d:")
                for i in range(len(d_out)):
                    if i < n:
                        if i == 0:
                            print("i=" + str(i+1) + ":  ", end=" ")
                            print(str(d_out[i]), end= ", ")
                        else:
                            if (i+1)%n == 0: 
                                print(str(d_out[i]))
                            else:
                                print(str(d_out[i]), end= ", ")
                    else:
                        if i % n == 0:
                            print("i=" + str((i)/n+1) + ":  ", end= " ")
                            print(str(d_out[i]), end= ", ")
                        else:
                            if (i+1)%n == 0: 
                                print(str(d_out[i]))    
                            else:
                                print(str(d_out[i]), end= ", ")
                        
                print()
                print("b:")
                for i in range(len(b_out)):
                    if i < n:
                        if i == 0:
                            print("i=" + str(i+1) + ":  ", end=" ")
                            print(str(b_out[i]), end= ", ")
                        else:
                            if (i+1)%n == 0: 
                                print(str(b_out[i]))
                            else:
                                print(str(b_out[i]), end= ", ")
                    else:
                        if i % n == 0:
                            print("i=" + str((i)/n+1) + ":  ", end= " ")
                            print(str(b_out[i]), end= ", ")
                        else:
                            if (i+1)%n == 0: 
                                print(str(b_out[i]))    
                            else:
                                print(str(b_out[i]), end= ", ")
    
                print()
                print("a0:")
                for i in range(n):
                    for j in range(n):
                        print(a0[i][j], end=", ")
                        if (j+1)%n == 0:
                            print()

                
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def main(n_in, m_in, input_file_in, I_in, n_ff_input):
    global n, m, start_time, I_input, nff, T_input
    start_time = time.time()
    input_file  = input_file_in
    n = n_in
    m = m_in
    I_input = I_in
    nff = n_ff_input
    createGraph(input_file)
    T_input = 20
    run()
    
if __name__ == "__main__":
    global instance
    folder_dataset = 'C:/Users/jgd/Documents/FireFighter/'
    dataset = [
        #['test.mtx' , 10, 9, [0], 1], # File name, n, m, fire sources, number of firefighters
        #['test2.mtx', 20, 19, [0,19], 1]
        ['grid4x4.mtx', 16, 24, [0], 1]
        #['grid3x3.mtx', 9, 12, [0], 1]
        #['test2.mtx', 20, 19, [0], 4]
        #['ba300_1.mtx' , 300, 299]
        #['ba100_1.mtx' , 100, 99]
        #['polbooks.mtx' , 105, 441]
        #['reds200_10.mtx' , 200, 380]
        #['lattice2d.mtx', 1089, 2112, [0], 1]
        #['dolphins.mtx', 62, 159, [0], 1]
        #['karate.mtx', 34, 78, [0], 1]
        #['lattice10x10.mtx', 100, 180, [55], 1]
        ]
    for i in range(len(dataset)):
        print("--------------------------------------------------------------")
        instance = dataset[i][0]
        print("instance: " + instance)
        main(dataset[i][1], dataset[i][2], folder_dataset + dataset[i][0], dataset[i][3], dataset[i][4])
        
