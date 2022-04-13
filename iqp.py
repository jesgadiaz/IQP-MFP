from gurobipy import *
import math
import numpy as np
import time
import networkx as nx
import random
import copy

def createGraph(input_file):
    global G, n, m, a0, start_time, apsp_time, n_nodes, n_edges, conn_comp, av_degree, density, card_comp, delta_m, mov
        
    print("graph created")
    for i in range(0, n):
        G.add_edge(i, i)
        
    if nx.is_connected(G):
        print("is connected")
        sp = dict(nx.all_pairs_shortest_path_length(G))
        delta_m = []
        for i in range(0,n):
            list = []
            for j in range(0,n):
                list.append(sp[i][j])
            delta_m.append(list)
        #for i in range(0, n):
        #    for j in range(0, n):
        #        if delta_m[i][j] == float("inf"):
        #            delta_m[i][j] = (n+1);
        del sp
    else:
        print("is not connected")
        delta_m = []
        for i in range(0,n):
            list = []
            for j in range(0,n):
                list.append(float("inf"))
            delta_m.append(list)
        connected_components = nx.connected_components(G)
        print(nx.number_connected_components(G))
        card_comp = []
        for component in connected_components:
            print("cc")
            # create subgraph    
            G_sub = nx.Graph()
            for v in component:
                G_sub.add_node(v)
                for e in G.edges:
                    if e[0]==v or e[1]==v:
                        G_sub.add_edge(e[0],e[1])
            card_comp.append(len(nx.nodes(G_sub)))
            
            sp = dict(nx.all_pairs_shortest_path_length(G_sub))
            for item1 in sp.items():
                for item2 in item1[1].items():
                    delta_m[item1[0]][item2[0]] = item2[1]
                    delta_m[item2[0]][item1[0]] = item2[1]
            del G_sub
            del sp
            
        for i in range(0, n):
            for j in range(0, n):
                if delta_m[i][j] == float("inf"):
                    delta_m[i][j] = (n+1)
                                    
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

def run(T_input, T2, Td, Tb):
    global size_best_solution, G, best_solution, model_feasible, T, total_runtime, runtime, n, m, a0, feasible, best_sequence, I_input, delta_m, feasible, best_defense, best_defense_time
    try:
        m = Model("mip1")
        m.Params.outputFlag = 1  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 2); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        
        T  = T_input # mid
        M  = n + 1
        
        
        # ------------------------------INPUT --------------------------------
        I = I_input  # Fire sources
        b0 = []
        d0 = []
        for i in range(n):
            if i in I:
                b0.append(1) # These are the vertices burned at time j=0
            else:
                b0.append(0)
        for i in range(n):
            if i == 0:
                d0.append(1) # The bulldozer begins at vertex 0
            else:
                d0.append(0)
                
        a = [] # Graph's adjacency matrix
        for i in range(n):
            list = []
            for j in range(n):
                list.append(0)
            a.append(list)
        for i in range(n):
            for j in range(n):
                if i in G.neighbors(j):
                    a[i][j] = 1
                    a[j][i] = 1
        # Matriz de movilidad del bulldozer
        #mov = []
        #for i in range(n):
        #    temp = []
        #    for j in range(n):
        #        temp.append(1.5)
        #    mov.append(temp)
            
        # Distancia de cada nodo hacia la fuente de fuego inicial
        dist_fuego = []
        for i in range(n):
            min_dist = float("inf")
            for j in I:
                if delta_m[i][j] < min_dist:
                    min_dist = delta_m[i][j]
            dist_fuego.append(min_dist)
            
        
        #---------------------------- VARIABLES -----------------------------------------------------------
        
        b = []
        for i in range(n):
            temp = []
            for j in range(Tb):
                temp.append(0)
            b.append(temp)
        for i in range(n):
            for j in range(Tb):
                b[i][j] = m.addVar(vtype=GRB.BINARY, name="b,%s" % str(i+1) + "," + str(j+1))
        
        d = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            d.append(temp)
        for i in range(n):
            for j in range(T):
                d[i][j] = m.addVar(vtype=GRB.BINARY, name="d,%s" % str(i+1) + "," + str(j+1))
                
        p = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            p.append(temp)
        for i in range(n):
            for j in range(T):
                p[i][j] = m.addVar(vtype=GRB.BINARY, name="p,%s" % str(i+1) + "," + str(j+1))
                

#---------------------------------------------------------------------------------                
        # Añado otra variable auxiliar p_prime
#        p_prime = []
#        for i in range(n):
#            temp = []
#            for j in range(T):
#                temp.append(0)
#            p_prime.append(temp)
#        for i in range(n):
#            for j in range(T):
#                p_prime[i][j] = m.addVar(vtype=GRB.BINARY, name="p_prime,%s" % str(i+1) + "," + str(j+1))                
#                
#        # Con esta garantizo que p_prime copie lo de p (los 1s)
#        for j in range(T):
#            for i in range(n):
#                m.addConstr(p_prime[i][j] >= p[i][j])
#                
#        # esta garantiza que haya sume a 1 las columnas de p_prime
#        for j in range(T):
#            sum_ = 0
#            for i in range(n):
#                sum_ = sum_ + p_prime[i][j]
#            m.addConstr(sum_ == 1)
#            
#        # Con esta, las columnas que suman cero en p, se convierten en p_prime = p_prime anterior
#        for j in range(T):
#            for i in range(n):
#                sum_ = 0
#                for k in range(n):
#                    sum_ = sum_ + p[k][j]
#                if j == 0:
#                    m.addConstr(p_prime[i][j] >= 0) # Deshabilitada, porque siempre defiende a alguien al inicio
#                else:
#                    m.addConstr(p_prime[i][j] >= p_prime[i][j-1] * (1 - sum_))
#---------------------------------------------------------------------------------                

                
        t = []
        for j in range(T):
            t.append(0)
        for j in range(T):
            t[j] = m.addVar(vtype=GRB.CONTINUOUS, name="t,%s" % str(j+1))
            
        t_floor = []
        for j in range(T):
            t_floor.append(0)
        for j in range(T):
            t_floor[j] = m.addVar(vtype=GRB.INTEGER, name="t_floor,%s" % str(j+1))
            
        b_prime = []
        for j in range(T):
            temp_1 = []
            for i in range(n):
                temp_2 = []
                for l in range(T2):
                    temp_2.append(0)
                temp_1.append(temp_2)
            b_prime.append(temp_1)
        for j in range(T):
            for i in range(n):
                for l in range(T2):
                    b_prime[j][i][l] = m.addVar(vtype=GRB.BINARY, name="b_prime,%s" % str(j+1) + "," + str(i+1) + "," + str(l+1))

        c = []
        for j in range(T):
            temp = []
            for i in range(T2):
                temp.append(0)
            c.append(temp)
        for j in range(T):
            for i in range(T2):
                c[j][i] = m.addVar(vtype=GRB.BINARY, name="c,%s" % str(j+1) + "," + str(i+1))
                    
        #---------------------------- CONSTRAINTS ---------------------------------------------------------
        
        for i in range(n): #---------------------------(2)
            for j in range(Tb):
                if j == 0:
                    m.addConstr(b[i][j] >= b0[i])
                else:
                    m.addConstr(b[i][j] >= b[i][j-1])
        
        for i in range(n): #-------------------------(3)
            for j in range(T):
                if j == 0:
                    m.addConstr(d[i][j] >= d0[i])
                else:
                    m.addConstr(d[i][j] >= d[i][j-1])

        for i in range(n): #---------------------------(4)
            for j in range(T):
                m.addConstr(b[i][j] + d[i][j] <= 1)
       
        
        for j in range(T): # Como conocemos T, se defienden a fuerzas
            sum_ = 0
            for i in range(n):
                sum_ = sum_ + p[i][j]
            m.addConstr(sum_ == 1)
     
#        for j in range(T): #---------------------------(5)
#            if j == 0:
#               sum_ = 0
#                for i in range(n):
#                    sum_ = sum_ + d[i][j] - d0[i]
#                m.addConstr(sum_ <= 1)
#            else:
#                sum_ = 0
#                for i in range(n):
#                    sum_ = sum_ + d[i][j] - d[i][j-1]
#                m.addConstr(sum_ <= 1)
            
#        sum_ = 0 # Esto es para obligar a que se defiendan nodos
#        for i in range(n):
#            sum_ = sum_ + b[i][T-1] + d[i][T-1]
#        m.addConstr(sum_ == n)

#        for j in range(1,T):
#            sum_1 = 0
#            sum_2 = 0
#            for i in range(n):
#                sum_1 = sum_1 + p[i][j]
#                sum_2 = sum_2 + p[i][j-1]
#            m.addConstr(sum_1 <= sum_2)
        

#        y = []
#        for i in range(n):
#            temp = []
#            for j in range(T):
#                temp.append(0)
#            y.append(temp)
#        for i in range(n):
#            for j in range(T):
#                y[i][j] = m.addVar(vtype=GRB.BINARY, name="y,%s" % str(i+1) + "," + str(j+1))#

#        z = []
#        for j in range(T):
#            z.append(0)
#        for j in range(T):
#            z[j] = m.addVar(vtype=GRB.BINARY, name="z,%s" % str(j+1))
            

#        # Hasta que no sea más posible, defiende un nuevo nodo
#        for j in range(T): #---------------------------(5 bis)
#            if j == 0:
#                sum_  = 0
#                sum_y = 0
#                for i in range(n):
#                    sum_ = sum_ + d[i][j] - d0[i]
#                    sum_y = sum_y + y[i][j]
#                    m.addConstr(y[i][j] >= d[i][j])            #OR
#                    m.addConstr(y[i][j] >= b[i][j])            #OR
#                    m.addConstr(y[i][j] <= d[i][j] + b[i][j])  #OR
#                    m.addConstr(z[j] <= y[i][j])            #AND
#                m.addConstr(z[j] >= sum_y - (n-1))      #AND
#                m.addConstr(sum_ >= 1 - z[j])
#                m.addConstr(sum_ <= (1 - z[j]) + z[j]*M)
#                #m.addConstr(sum_ == 1)
#            else:
#                sum_ = 0
#                sum_y = 0
#                for i in range(n):
#                    sum_ = sum_ + d[i][j] - d[i][j-1]
#                    sum_y = sum_y + y[i][j]
#                    m.addConstr(y[i][j] >= d[i][j])            #OR
#                    m.addConstr(y[i][j] >= b[i][j])            #OR
#                    m.addConstr(y[i][j] <= d[i][j] + b[i][j])  #OR
#                    m.addConstr(z[j] <= y[i][j])            #AND
#                m.addConstr(z[j] >= sum_y - (n-1))      #AND
#                m.addConstr(sum_ >= 1 - z[j])
#                m.addConstr(sum_ <= (1 - z[j]) + z[j]*M)
                #m.addConstr(sum_ == 1)
                
                
                

        for i in range(n): #-------------------------(6)
            for j in range(T):
                if j == 0:
                    m.addConstr(p[i][j] == d[i][j] - d0[i])
                else:
                    m.addConstr(p[i][j] == d[i][j] - d[i][j-1])
        
                    
        for j in range(T): #-----------------------(7) 
            dist_nodos = 0
            for i in range(n):
                sum_1 = 0
                if j == 0:
                    for k in range(n):
                        sum_1 = sum_1 + p[k][j] * mov[k][i]
                    dist_nodos = dist_nodos + sum_1 * d0[i]
                else:
                    for k in range(n):
                        sum_1 = sum_1 + p[k][j] * mov[k][i]
                    dist_nodos = dist_nodos + sum_1 * p[i][j-1]
            if j == 0:
                m.addConstr(t[j] == dist_nodos)
            else:
                m.addConstr(t[j] == dist_nodos + t[j-1])
                    
        for j in range(1,T): #----------------------------(8 y 9)
            for i in range(n):
                for k in range(T2):
                    if k == 0:
                        m.addConstr(b_prime[j][i][k] >= b[i][j-1])
                    else:
                        m.addConstr(b_prime[j][i][k] >= b_prime[j][i][k-1])

        for j in [0]: #----------------------------(8 y 9) con j=0
            for i in range(n):
                for k in range(T2):
                    if k == 0:
                        m.addConstr(b_prime[j][i][k] >= b0[i])
                    else:
                        m.addConstr(b_prime[j][i][k] >= b_prime[j][i][k-1])
 
                        
 
 ### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OJO ACÁ !!!!!!!!!!!!!!
    
        # Se debe defender alguien desde el inicio
#        sum_ = 0
#        for i in range(n):
#            sum_ = sum_ + p[i][0]
#        m.addConstr(sum_ == 1)
        # Y se deben defender siempre en cada ronda (consecutivamente)
#        for j in range(1,T):
#            sum_1 = 0
#            sum_2 = 0
#            for i in range(n):
#                sum_1 = sum_1 + p[i][j-1]
#                sum_2 = sum_2 + p[i][j]
#            m.addConstr(sum_2 <= sum_1)
 
        # Manejar T como entrada del problema y checar la solución
        # para ver si es consistente y minimiza los quemados
 
        #for j in range(T):
        for j in range(T): # Hay que encontrar este valor en automático
            sum_1 = 0
            sum_2 = 0
            for i in range(n):
                sum_1 = sum_1 + p[i][j]
                sum_2 = sum_2 + dist_fuego[i] * p[i][j]
            m.addConstr(sum_2 >= t[j] * sum_1)
 
### !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! OJO ACÁ !!!!!!!!!!!!!!    
 
    
 
    
 
    
        # w[i][j]=1 si el nodo i se puede defender en el tiempo j
#        w = []
#        for i in range(n):
#            temp = []
#            for j in range(T):
#                temp.append(0)
#            w.append(temp)
#        for i in range(n):
#            for j in range(T):
#                w[i][j] = m.addVar(vtype=GRB.BINARY, name="w,%s" % str(i+1) + "," +str(j+1))
#        for i in range(n):
#            for j in range(T):
#                # w[i][j] = 1 si su distancia es mayor o igual a t[j]
#                m.addConstr(dist_fuego[i] >= t[j] - M * (1-w[i][j]))
#                m.addConstr(t[j] >= dist_fuego[i] + 0.0001 - M * w[i][j]) #debe ser 0.0001
                
        # Se prohibe que se quemen los nodos i en el tiempo j
        # Esto debe desactivarse después del tiempo max_t, el cual es la distancia 
        # desde la fuente de quema hasta el nodo más alejado diferente del almacén
#        max_dist = 0
#        for i in range(1,n):
#            min_dist = float("inf")
#            for j in I:
#                if delta_m[i][j] < min_dist:
#                    min_dist = delta_m[i][j]
#            if min_dist > max_dist:
#                max_dist = min_dist
#        for j in range(max_dist): # Es hasta la columna con t[j] >= max_dist
#            for i in range(n):    
#                m.addConstr(p[i][j] <= w[i][j])
                
                                
        for j in range(1,T): #----------------------------(10 y 11)
            for i in range(n):
                for l in range(T2):
                    if l == 0:
                        for k in G.neighbors(i):
                            m.addConstr(b_prime[j][i][l] + d[i][j] >= b[k][j-1] * c[j][l])
                    else:
                        for k in G.neighbors(i):
                            m.addConstr(b_prime[j][i][l] + d[i][j] >= b_prime[j][k][l-1] * c[j][l])

        for j in [0]: #----------------------------(10 y 11) con j=0
            for i in range(n):
                for l in range(T2):
                    if l == 0:
                        for k in G.neighbors(i):
                            m.addConstr(b_prime[j][i][l] + d[i][j] >= b0[k] * c[j][l])
                    else:
                        for k in G.neighbors(i):
                            m.addConstr(b_prime[j][i][l] + d[i][j] >= b_prime[j][k][l-1] * c[j][l])
                    
                            


#----- Acá se evita que se quemen nodos cuando c[j][i] = 0---------------------------------------------------------------------
        for j in range(1,T):
            for i in range(n):
                for l in range(T2):
                    if l == 0:
                        m.addConstr(b_prime[j][i][l] - b[i][j-1] <= c[j][l])
                    else:
                        m.addConstr(b_prime[j][i][l] - b_prime[j][i][l-1] <= c[j][l])

        for j in [0]:
            for i in range(n):
                for l in range(T2):
                    if l == 0:
                        m.addConstr(b_prime[j][i][l] - b0[i] <= c[j][l])
                    else:
                        m.addConstr(b_prime[j][i][l] - b_prime[j][i][l-1] <= c[j][l])
#--------------------------------------------------------------------------

                    
                    
#--------------------------------------------------------------------------------
# Esto es para evitar que los nodos que no tienen vecinos con fuego, se quemen
# Checar si voltear i y l
        for j in range(1,T): #----------------------------()
            for i in range(n):
                for l in range(T2):
                    sum_ = 0
                    if l == 0:
                        for k in G.neighbors(i):
                            sum_ = sum_ + b[k][j-1]
                    else:
                        for k in G.neighbors(i):
                            sum_ = sum_ + b_prime[j][k][l-1]
                    m.addConstr(b_prime[j][i][l] <= sum_)

        for j in [0]: #----------------------------() con j=0
            for i in range(n):
                for l in range(T2):
                    sum_ = 0
                    if l == 0:
                        for k in G.neighbors(i):
                            sum_ = sum_ + b0[k]
                    else:
                        for k in G.neighbors(i):
                            sum_ = sum_ + b_prime[j][k][l-1]
                    m.addConstr(b_prime[j][i][l] <= sum_)
#--------------------------------------------------------------------------------                            
                    

                    
#        for j in range(1,T): #----------------------------(12 y 13)
#            for i in range(T2):
#                sum_ = 0
#                for k in range(n):
#                    if j == 0:
#                        m.addConstr(c[j][i] >= b_prime[j][k][i] - b[k][j-1])
#                        sum_ = sum_ + b_prime[j][k][i] - b[k][j-1]
#                    else:
#                        m.addConstr(c[j][i] >= b_prime[j][k][i] - b_prime[j][k][i-1])
#                        sum_ = sum_ + b_prime[j][k][i] - b_prime[j][k][i-1]
#                m.addConstr(c[j][i] <= sum_)
                
        #----------------------------(12 y 13) con j=0 ????
#        for i in range(T2):
#            sum_ = 0
#            for k in range(n):
#                if i == 0:
#                    m.addConstr(c[0][i] >= b_prime[0][k][i] - b0[k])
#                    sum_ = sum_ + b_prime[0][k][i] - b0[k]
#                else:
#                    m.addConstr(c[0][i] >= b_prime[0][k][i] - b_prime[0][k][i-1])
#                    sum_ = sum_ + b_prime[0][k][i] - b_prime[0][k][i-1]
#            m.addConstr(c[0][i] <= sum_)



                
        for j in range(T): #--------------------------(14)  ????
            sum_ = 0
            for i in range(T2):
                sum_ = sum_ + c[j][i]
            if j == 0:
                m.addConstr(sum_ == t_floor[j])
            else:
                m.addConstr(sum_ == t_floor[j] - t_floor[j-1])
                
        for j in range(T): # ----------------------(14 bis)
            for i in range(1,T2):
                m.addConstr(c[j][i] <= c[j][i-1])
                    
        epsilon = 0.001
        for j in range(T):
            m.addConstr(t_floor[j] <= t[j])
            m.addConstr(t_floor[j] + 1 >= t[j] + epsilon)
            
        for j in range(T): # -----------------------(15)
            for i in range(n):
                m.addConstr(b[i][j] == b_prime[j][i][T2-1])
                
                
        #-------------------------------------------------
        # De T+1 hasta Tb se queman lo nodos vecinos de manera normal
        for j in range(T,Tb-1):
            for i in range(n):
                for k in G.neighbors(i):
                    m.addConstr(b[k][j] >= b[i][j-1] * (1-d[k][T-1]))
                    
                
        #m.addConstr(d[8][1] == 1)
                                        
        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------        
        b_transpose = np.array(b).T.tolist()
        m.setObjective(sum(b_transpose[Tb-1]), GRB.MINIMIZE)#-----------------------------------------------(1)
        #d_transpose = np.array(d).T.tolist()
        #m.setObjective(sum(d_transpose[T-1]), GRB.MAXIMIZE)#-----------------------------------------------(1)
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        #m.computeIIS()
        #m.write("model.lp")
        
        
        m.optimize()
        runtime = m.Runtime
        print("The run time is %f" % runtime)
        model_feasible = True
        feasible = True
        
        if m.status == GRB.INFEASIBLE:
            model_feasible = False
        else:
            
            print("Obj:", m.objVal)
            
            b_out       = []
            d_out       = []
            p_out       = []
            c_out       = []
            w_out       = []
            t_out       = []
            t_floor_out = []
            b_prime_out = []
            p_prime_out = []
            for v in m.getVars():
                varName = v.varName
                varNameSplit = varName.split(',')
                if varNameSplit[0] == 'd':
                   d_out.append(v.x)
                if varNameSplit[0] == 'b':
                   b_out.append(v.x)
                if varNameSplit[0] == 'b_prime':
                   b_prime_out.append(v.x)
                if varNameSplit[0] == 'p_prime':
                   p_prime_out.append(v.x)
                if varNameSplit[0] == 't':
                   t_out.append(v.x)
                if varNameSplit[0] == 't_floor':
                   t_floor_out.append(v.x)
                if varNameSplit[0] == 'p':
                   p_out.append(v.x)
                if varNameSplit[0] == 'c':
                   c_out.append(v.x)
                if varNameSplit[0] == 'w':
                   w_out.append(v.x)          
            
            t = []
            for i in range(len(t_out)):
                t.append(t_out[i])
    
            t_floor = []
            for i in range(len(t_floor_out)):
                t_floor.append(t_floor_out[i])
    
            c = []
            for i in range(len(c_out)):
                c.append(c_out[i])
                   
            d = []
            for i in range(n):
                temp = []
                for j in range(T):
                    temp.append(0)
                d.append(temp)                
            j = 0
            k = 0
            
            for i in range(len(d_out)):
                if d_out[i] > 0.9:
                    d[j][k] = 1
                else:
                    d[j][k] = 0
                #d[j][k] = d_out[i]
                k =  k + 1
                if (i+1) % T == 0:
                    j = j + 1
                    k = 0
    
    
            w = []
            for i in range(n):
                temp = []
                for j in range(T):
                    temp.append(0)
                w.append(temp)                
            j = 0
            k = 0
            for i in range(len(w_out)):
                if w_out[i] > 0.9:
                    w[j][k] = 1
                else:
                    w[j][k] = 0
                #w[j][k] = w_out[i]
                k =  k + 1
                if (i+1) % T == 0:
                    j = j + 1
                    k = 0
    
    
            p = []
            for i in range(n):
                temp = []
                for j in range(T):
                    temp.append(0)
                p.append(temp)                
            j = 0
            k = 0
            
            for i in range(len(p_out)):
                if p_out[i] > 0.9:
                    p[j][k] = 1
                else:
                    p[j][k] = 0
                #p[j][k] = p_out[i]
                k =  k + 1
                if (i+1) % T == 0:
                    j = j + 1
                    k = 0
    
    
            p_prime = []
            for i in range(n):
                temp = []
                for j in range(T):
                    temp.append(0)
                p_prime.append(temp)                
            j = 0
            k = 0
            
            for i in range(len(p_prime_out)):
                if p_prime_out[i] > 0.9:
                    p_prime[j][k] = 1
                else:
                    p_prime[j][k] = 0
                #p[j][k] = p_out[i]
                k =  k + 1
                if (i+1) % T == 0:
                    j = j + 1
                    k = 0
    
    
            b = []
            for i in range(n):
                temp = []
                for j in range(Tb):
                    temp.append(0)
                b.append(temp)
            j = 0
            k = 0
            
            for i in range(len(b_out)):
                if b_out[i] > 0.9:
                    b[j][k] = 1
                else:
                    b[j][k] = 0
                #b[j][k] = b_out[i]
                k =  k + 1
                if (i+1) % Tb == 0:
                    j = j + 1
                    k = 0
                    
            b_prime = []
            for j in range(T):
                temp1 = []
                for i in range(n):
                    temp2 = []
                    for j in range(T2):
                        temp2.append(0)
                    temp1.append(temp2)
                b_prime.append(temp1)
    
            j = 0
            k = 0
            l = 0
            for i in range(len(b_prime_out)):
                if b_prime_out[i] > 0.9:
                    b_prime[l][j][k] = 1
                else:
                    b_prime[l][j][k] = 0
                k =  k + 1
                if (i+1) % T2 == 0:
                    j = j + 1
                    k = 0  
                if (i+1) % (n*T2) == 0:
                    j = 0
                    k = 0
                    l += 1
                    
            # Defense plan
            defense = []
            time    = []
            for j in range(T):
                k = -1
                for i in range(n):
                   if p[i][j] == 1:
                       k = i
                if k != -1:
                    defense.append(k)
                    time.append(t[j])
            defense_out = copy.deepcopy(defense)
            time_out = copy.deepcopy(time)
            #print("Defense plan: " + str(defense))
            #print("Defense plan: " + str(time))
            # Check if the solution is feasible
            #feasible = True
            # Fuentes iniciales de fuego
            burned = [0 for i in range(n)]
            for i in I:
                burned[i] = 1
            next_burned = copy.deepcopy(burned)
            # Radio de I
            max_dist = 0
            for i in range(1,n):
                min_dist = float("inf")
                for j in I:
                    if delta_m[i][j] < min_dist:
                        min_dist = delta_m[i][j]
                if min_dist > max_dist:
                    max_dist = min_dist            
            # Empezamos con las rondas de quema
            for j in range(1,max_dist+1):
                #print(burned)
                for i in range(n):
                    if burned[i] == 1:
                        for k in G.neighbors(i):
                            if k in defense:
                                #print(burned)
                                #print(defense)
                                #print(time)
                                index = defense.index(k)
                                if time[index] > j:
                                    feasible = False
                                #sacar de defense porque ya lo checamos
                                index = defense.index(k)
                                defense.pop(index)
                                time.pop(index)
                            else:
                                if k not in defense_out:
                                    next_burned[k] = 1
                burned = copy.deepcopy(next_burned)
            
            if feasible:
                if sum(burned) < size_best_solution:
                    size_best_solution = sum(burned)
                    best_solution = copy.deepcopy([defense_out, time_out])
                    
            #if feasible:
            #    best_defense      = defense_out[:Td]
            #    best_defense_time = time_out[:Td]
            #print(best_defense)
            #print(best_defense_time)
                
                    
            #print("d ---------------------------")
            #for i in d:
            #    print(i)
            #print("p ---------------------------")
            #for i in p:
            #    print(i)
            #print("p_prime ---------------------------")
            #for i in p_prime:
            #    print(i)
            #print("d0 ---------------------------")            
            #print(d0)
            #print("b ---------------------------")
            #for i in b:
            #    print(i)
            #print("w ---------------------------")
            #for i in w:
            #    print(i)            
            #print("b0 ---------------------------")            
            #print(b0)
            #print("t ---------------------------")            
            #print(t)
            #print("t_floor ---------------------")            
            #print(t_floor) 
            #v = 0
            #for i in range(len(b_prime)):
            #for i in [0,1]:
            #print()
            #    print("b_prime_" + str(i))
            #    for j in b_prime[i]:
            #        print(j)
            #    for k in range(T2):
            #        print(str(int(c[v])),end=", ")
            #        v+=1
    except GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

def main(n_in, m_in, input_file_in, I_in):
    global n, m, start_time, I_input#, T_input, T2
    start_time = time.time()
    input_file  = input_file_in
    n = n_in
    m = m_in
    I_input = I_in
    createGraph(input_file)
    #T_input = 14
    #T2 = 5
    
if __name__ == "__main__":
    global size_best_solution, best_solution, model_feasible, instance, feasible, best_defense, best_defense_time, runtime, G, mov
    folder_dataset = 'C:/Users/jgd/Documents/FireFighter/'
    dataset = [
        #['lattice4x4.mtx', 16, 24, [15]]
        #['AdjList_1.adjlist' ,[18]]
        #['AdjList_1_new.adjlist' ,[26]]
        ['AdjList_2_new.adjlist' ,[28]]
        #['AdjList_3_new.adjlist' ,[17]]
        #['Tree_new.adjlist' ,[18]]
        #['lattice2x2.mtx', 4, 4, [3]]
        #['path4.mtx', 4, 3, [3]]
        #['test_lattice2x2.mtx', 5, 5, [2]]
        #['test2_lattice2x2.mtx', 6, 6, [2]]
        #['test3_lattice2x2.mtx', 8, 8, [2]]
        #['path6.mtx', 6, 5, [5]]
        ]
    # Tree
    #mov = [[0.0, 2.69453808, 3.64161647, 12.34325843, 2.86009481, 6.1374802, 5.32246651, 11.29659716, 3.31525444, 0.78519918, 12.41331877, 12.71590282, 10.85152542, 5.30531019, 14.62988656, 15.66729075, 5.1738933, 3.58142053, 8.10023862, 7.72959113, 2.08169749, 4.59862711, 6.84565681, 6.91266051, 16.3123257, 4.8405079], [2.69453808, 0.0, 6.33286266, 9.64875876, 4.29425366, 3.62836634, 6.02264541, 8.60592167, 3.06480263, 2.25476522, 9.73701553, 10.12917852, 8.17368452, 7.60734617, 12.13698276, 13.27478887, 2.47957666, 5.94352464, 5.43828894, 5.04438059, 4.76417279, 2.78090837, 7.02480049, 7.70310563, 13.81340365, 6.77107657], [3.64161647, 6.33286266, 0.0, 15.9780014, 4.00804448, 9.64134273, 6.17605611, 14.93774382, 5.90190981, 4.16710899, 16.01513254, 16.21524974, 14.45748061, 3.5593164, 18.01879751, 18.94081739, 8.81211911, 2.56523183, 11.73700783, 11.34952072, 1.57468355, 7.8158263, 7.95886858, 7.21101714, 19.69790598, 4.25692634], [12.34325843, 9.64875876, 15.9780014, 0.0, 13.18906052, 6.65734739, 13.53490248, 1.21019963, 11.06493431, 11.83677242, 1.21847207, 3.11832915, 1.83175951, 16.84611761, 5.18439858, 6.82614445, 7.17018836, 15.33762043, 4.44343542, 4.65384074, 14.40544845, 8.8601632, 13.28726193, 14.79876271, 6.2995567, 15.67399502], [2.86009481, 4.29425366, 4.00804448, 13.18906052, 0.0, 6.54302706, 2.55203149, 12.30111758, 2.2016778, 2.46256654, 12.99438637, 12.87478936, 11.49295696, 7.07799438, 14.46175496, 15.24009862, 6.40054883, 5.49462281, 9.40801035, 8.57909014, 2.92421481, 4.40130551, 4.23596374, 4.06829594, 16.11947556, 7.11223997], [6.1374802, 3.62836634, 9.64134273, 6.65734739, 6.54302706, 0.0, 7.02023245, 5.83650559, 4.40798494, 5.49154145, 6.46790425, 6.58377272, 4.9515114, 11.23322217, 8.51841325, 9.64893718, 1.89886861, 9.55770251, 3.47662582, 2.17386813, 8.07913676, 2.21692624, 7.13140097, 8.44665983, 10.19802227, 10.37950982], [5.32246651, 6.02264541, 6.17605611, 13.53490248, 2.55203149, 7.02023245, 0.0, 12.82257679, 3.01623754, 4.77488592, 13.11344175, 12.62301977, 11.7338606, 9.50479511, 13.85853526, 14.38249205, 7.54307763, 7.98929774, 10.34049008, 9.1927582, 5.38282355, 4.84396878, 1.78281258, 1.68452975, 15.44610294, 9.64217133], [11.29659716, 8.60592167, 14.93774382, 1.21019963, 12.30111758, 5.83650559, 12.82257679, 0.0, 10.2247338, 10.82382343, 2.09441163, 3.84432566, 1.64695079, 15.6963289, 6.09367811, 7.72561094, 6.12954017, 14.21040127, 3.28797438, 3.72233099, 13.37009259, 8.05330675, 12.71380328, 14.16409044, 7.34540194, 14.50372216], [3.31525444, 3.06480263, 5.90190981, 11.06493431, 2.2016778, 4.40798494, 3.01623754, 10.2247338, 0.0, 2.55091774, 10.8229774, 10.67314007, 9.33868012, 8.46980964, 12.27957352, 13.1003455, 4.58428986, 6.76531945, 7.4908784, 6.51112833, 4.51501243, 2.21931093, 3.98174339, 4.67500764, 13.94327444, 8.14568881], [0.78519918, 2.25476522, 4.16710899, 11.83677242, 2.46256654, 5.49154145, 4.77488592, 10.82382343, 2.55091774, 0.0, 11.84912271, 12.07478086, 10.29278088, 6.0790813, 13.94477075, 14.95145432, 4.69860935, 4.35512136, 7.68864491, 7.1942034, 2.59516155, 3.84654483, 6.20912235, 6.41300434, 15.6272628, 5.62226734], [12.41331877, 9.73701553, 16.01513254, 1.21847207, 12.99438637, 6.46790425, 13.11344175, 2.09441163, 10.8229774, 11.84912271, 0.0, 1.90003438, 1.56357422, 17.13358255, 4.02199744, 5.66167899, 7.29470888, 15.5678518, 4.89030374, 4.69686133, 14.4410741, 8.60373002, 12.71754856, 14.28325341, 5.25445843, 16.03227189], [12.71590282, 10.12917852, 16.21524974, 3.11832915, 12.87478936, 6.58377272, 12.62301977, 3.84432566, 10.67314007, 12.07478086, 1.90003438, 0.0, 2.6118074, 17.70147714, 2.3518661, 3.93205427, 7.83263402, 16.06600072, 5.99873107, 5.31296985, 14.65862467, 8.50059261, 11.98375994, 13.61875396, 3.86709812, 16.71557835], [10.85152542, 8.17368452, 14.45748061, 1.83175951, 11.49295696, 4.9515114, 11.7338606, 1.64695079, 9.33868012, 10.29278088, 1.56357422, 2.6118074, 0.0, 15.58271487, 4.96315262, 6.53297558, 5.73233014, 14.00948474, 3.45819553, 3.13352906, 12.88307769, 7.12195452, 11.45710668, 12.9758221, 6.42125175, 14.49783156], [5.30531019, 7.60734617, 3.5593164, 16.84611761, 7.07799438, 11.23322217, 9.50479511, 15.6963289, 8.46980964, 6.0790813, 17.13358255, 17.70147714, 15.58271487, 0.0, 19.74369024, 20.86788685, 9.8867652, 1.72409828, 12.40949609, 12.47385352, 4.15383123, 9.89548123, 11.26718792, 10.69383476, 21.41887369, 1.50735676], [14.62988656, 12.13698276, 18.01879751, 5.18439858, 14.46175496, 8.51841325, 13.85853526, 6.09367811, 12.27957352, 13.94477075, 4.02199744, 2.3518661, 4.96315262, 19.74369024, 0.0, 1.64194246, 9.96019735, 18.07538198, 8.33218791, 7.53462099, 16.49043977, 10.20833345, 12.95956846, 14.63838209, 1.68260264, 18.82844936], [15.66729075, 13.27478887, 18.94081739, 6.82614445, 15.24009862, 9.64893718, 14.38249205, 7.72561094, 13.1003455, 14.95145432, 5.66167899, 3.93205427, 6.53297558, 20.86788685, 1.64194246, 0.0, 11.22266398, 19.17578297, 9.82799986, 8.90862517, 17.44571855, 11.13886346, 13.31099511, 14.99462349, 1.51257476, 20.02425449], [5.1738933, 2.47957666, 8.81211911, 7.17018836, 6.40054883, 1.89886861, 7.54307763, 6.12954017, 4.58428986, 4.69860935, 7.29470888, 7.83263402, 5.73233014, 9.8867652, 9.96019735, 11.22266398, 0.0, 8.28568492, 3.02503932, 2.60056235, 7.24217467, 2.90742183, 8.06758082, 9.13380464, 11.61420387, 8.88477454], [3.58142053, 5.94352464, 2.56523183, 15.33762043, 5.49462281, 9.55770251, 7.98929774, 14.21040127, 6.76531945, 4.35512136, 15.5678518, 16.06600072, 14.00948474, 1.72409828, 18.07538198, 19.17578297, 8.28568492, 0.0, 10.92506501, 10.88435265, 2.61067234, 8.17412445, 9.72321638, 9.28935165, 19.75367659, 1.74174469], [8.10023862, 5.43828894, 11.73700783, 4.44343542, 9.40801035, 3.47662582, 10.34049008, 3.28797438, 7.4908784, 7.68864491, 4.89030374, 5.99873107, 3.45819553, 12.40949609, 8.33218791, 9.82799986, 3.02503932, 10.92506501, 0.0, 1.58089135, 10.18172332, 5.50726129, 10.59249795, 11.84585271, 9.85510696, 11.23092156], [7.72959113, 5.04438059, 11.34952072, 4.65384074, 8.57909014, 2.17386813, 9.1927582, 3.72233099, 6.51112833, 7.1942034, 4.69686133, 5.31296985, 3.13352906, 12.47385352, 7.53462099, 8.90862517, 2.60056235, 10.88435265, 1.58089135, 0.0, 9.77500048, 4.369047, 9.26792185, 10.61550814, 9.1523327, 11.43006949], [2.08169749, 4.76417279, 1.57468355, 14.40544845, 2.92421481, 8.07913676, 5.38282355, 13.37009259, 4.51501243, 2.59516155, 14.4410741, 14.65862467, 12.88307769, 4.15383123, 16.49043977, 17.44571855, 7.24217467, 2.61067234, 10.18172332, 9.77500048, 0.0, 6.3072182, 7.1281481, 6.68967693, 18.17146896, 4.30327841], [4.59862711, 2.78090837, 7.8158263, 8.8601632, 4.40130551, 2.21692624, 4.84396878, 8.05330675, 2.21931093, 3.84654483, 8.60373002, 8.50059261, 7.12195452, 9.89548123, 10.20833345, 11.13886346, 2.90742183, 8.17412445, 5.50726129, 4.369047, 6.3072182, 0.0, 5.17553967, 6.34196185, 11.88525069, 9.28096237], [6.84565681, 7.02480049, 7.95886858, 13.28726193, 4.23596374, 7.13140097, 1.78281258, 12.71380328, 3.98174339, 6.20912235, 12.71754856, 11.98375994, 11.45710668, 11.26718792, 12.95956846, 13.31099511, 8.06758082, 9.72321638, 10.59249795, 9.26792185, 7.1281481, 5.17553967, 0.0, 1.68452525, 14.47827637, 11.3473849], [6.91266051, 7.70310563, 7.21101714, 14.79876271, 4.06829594, 8.44665983, 1.68452975, 14.16409044, 4.67500764, 6.41300434, 14.28325341, 13.61875396, 12.9758221, 10.69383476, 14.63838209, 14.99462349, 9.13380464, 9.28935165, 11.84585271, 10.61550814, 6.68967693, 6.34196185, 1.68452525, 0.0, 16.16189478, 10.99268415], [16.3123257, 13.81340365, 19.69790598, 6.2995567, 16.11947556, 10.19802227, 15.44610294, 7.34540194, 13.94327444, 15.6272628, 5.25445843, 3.86709812, 6.42125175, 21.41887369, 1.68260264, 1.51257476, 11.61420387, 19.75367659, 9.85510696, 9.1523327, 18.17146896, 11.88525069, 14.47827637, 16.16189478, 0.0, 20.49131732], [4.8405079, 6.77107657, 4.25692634, 15.67399502, 7.11223997, 10.37950982, 9.64217133, 14.50372216, 8.14568881, 5.62226734, 16.03227189, 16.71557835, 14.49783156, 1.50735676, 18.82844936, 20.02425449, 8.88477454, 1.74174469, 11.23092156, 11.43006949, 4.30327841, 9.28096237, 11.3473849, 10.99268415, 20.49131732, 0.0]]
    # AdjList_1
    #mov = [0.0, 13.87300946, 13.03937355, 10.26931471, 8.88709628, 7.15949726, 15.79267074, 20.34855371, 9.32694321, 8.11218109, 9.593138, 5.55936489, 2.69133386, 12.64895326, 11.83495959, 11.70879828, 6.4373446, 7.78414476, 6.68451163, 7.79300168, 17.34996024, 17.60395722, 9.79397585, 9.85718622, 19.21619553, 3.82862739, 8.7164367, 13.18288407, 17.53956817, 11.70369435, 8.67701736] ,[13.87300946, 0.0, 4.24323817, 3.67236146, 5.7303412, 8.03122229, 2.23783709, 6.68986255, 8.04866435, 7.87049496, 14.98945227, 10.18991545, 12.93697881, 5.95639905, 2.40239869, 6.07824264, 11.11678063, 8.44551052, 8.60010133, 11.75839269, 6.05118257, 5.30404834, 6.63675716, 5.70534735, 5.72297753, 11.6685413, 13.59625527, 5.22122152, 4.31369829, 4.7267962, 10.82099606] ,[13.03937355, 4.24323817, 0.0, 5.07728371, 7.30347763, 9.24700528, 6.18408462, 9.87588616, 4.83585554, 5.33980017, 16.82365268, 11.04368349, 12.94844797, 1.75205743, 2.88335684, 2.16997718, 12.37288352, 10.07524752, 6.58401409, 13.4353958, 10.28382056, 9.44039185, 3.62189896, 7.99620437, 9.22408933, 11.99509716, 15.37513181, 0.98681961, 8.18963269, 1.3594828, 12.82301273] ,[10.26931471, 3.67236146, 5.07728371, 0.0, 2.323665, 4.53116383, 5.52837509, 10.09646867, 6.33673621, 5.56947424, 11.91344487, 6.61560053, 9.26589665, 6.15052067, 2.34510263, 5.70239291, 7.69582198, 5.1626248, 5.78125939, 8.55787492, 7.88908566, 7.73356999, 5.10321362, 2.92762016, 8.99537615, 8.01345658, 10.47528224, 5.89713902, 7.36247894, 4.48740453, 7.82270527] ,[8.88709628, 5.7303412, 7.30347763, 2.323665, 0.0, 2.30142214, 7.21729467, 11.64918299, 7.50071967, 6.39406928, 9.59251253, 4.47589495, 7.43633595, 8.17549729, 4.66479455, 7.58322465, 5.41598382, 2.83912446, 6.06284458, 6.23423946, 8.50188127, 8.71843188, 6.53196685, 1.30475623, 10.46645486, 6.07234356, 8.15276462, 8.05259831, 8.75921413, 6.53127564, 5.52542753] ,[7.15949726, 8.03122229, 9.24700528, 4.53116383, 2.30142214, 0.0, 9.47058496, 13.82802286, 8.39464982, 7.08258705, 7.59649398, 2.19285638, 5.33499761, 9.85930991, 6.81434145, 9.13965091, 3.17264835, 1.12205638, 6.29177839, 4.2081104, 10.28569975, 10.69209409, 7.73901619, 2.81956828, 12.62381157, 3.90893715, 6.14837579, 9.88932377, 10.90755626, 8.29159556, 3.82122721] ,[15.79267074, 2.23783709, 6.18408462, 5.52837509, 7.21729467, 9.47058496, 0.0, 4.57337959, 10.27306186, 10.10601659, 15.83151905, 11.66343562, 14.62875912, 7.93540704, 4.63742738, 8.17827032, 12.37857421, 9.66627221, 10.80194322, 12.78714438, 4.30212228, 3.27891144, 8.86209542, 6.82343432, 3.52381635, 13.28684489, 14.50755284, 7.11414623, 2.08216715, 6.86896002, 11.70959286] ,[20.34855371, 6.68986255, 9.87588616, 10.09646867, 11.64918299, 13.82802286, 4.57337959, 0.0, 14.45411313, 14.47312945, 19.51470416, 16.00751228, 19.08541924, 11.55480272, 9.05636221, 12.03251743, 16.55079062, 13.86878021, 15.28080044, 16.74021478, 5.10872591, 3.83728518, 13.06943212, 11.06218161, 1.24056048, 17.70810142, 18.29819095, 10.64149624, 2.93918273, 10.87722827, 15.55557199] ,[9.32694321, 8.04866435, 4.83585554, 6.33673621, 7.50071967, 8.39464982, 10.27306186, 14.45411313, 0.0, 1.40630666, 15.58996582, 9.37227373, 10.07329709, 3.69386905, 5.70986429, 2.81169309, 10.9383549, 9.48981659, 2.75163283, 12.34462784, 13.79041051, 13.27898131, 1.41196102, 8.70828004, 13.65033066, 9.5787741, 14.19733006, 4.48787087, 12.35512155, 3.61488034, 12.19699495] ,[8.11218109, 7.87049496, 5.33980017, 5.56947424, 6.39406928, 7.08258705, 10.10601659, 14.47312945, 1.40630666, 0.0, 14.19540853, 7.97688787, 8.70670776, 4.58778427, 5.46866104, 3.6463882, 9.54749941, 8.1910259, 1.43237716, 10.96807515, 13.31152777, 12.93395999, 1.75956695, 7.65194891, 13.57881789, 8.17796847, 12.80860634, 5.23515614, 12.16922105, 3.99564885, 10.86110176] ,[9.593138, 14.98945227, 16.82365268, 11.91344487, 9.59251253, 7.59649398, 15.83151905, 19.51470416, 15.58996582, 14.19540853, 0.0, 6.21852189, 6.99240744, 17.44387647, 14.25705759, 16.69623841, 4.65581362, 6.76328816, 13.06730703, 3.39012822, 14.82726971, 15.79334966, 15.16990059, 9.2924934, 18.2779851, 6.51126118, 1.44894598, 17.48502524, 16.64965172, 15.88775764, 4.16856321] ,[5.55936489, 10.18991545, 11.04368349, 6.61560053, 4.47589495, 2.19285638, 11.66343562, 16.00751228, 9.37227373, 7.97688787, 6.21852189, 0.0, 3.3255759, 11.42922451, 8.8175326, 10.62470016, 1.58412805, 2.34615336, 6.87329623, 3.09159006, 12.31382888, 12.80286736, 9.00424136, 4.95939238, 14.7979951, 1.87245486, 4.84219825, 11.58558897, 13.08101108, 9.96271172, 3.44220958] ,[2.69133386, 12.93697881, 12.94844797, 9.26589665, 7.43633595, 5.33499761, 14.62875912, 19.08541924, 10.07329709, 8.70670776, 6.99240744, 3.3255759, 0.0, 12.92102814, 11.20494726, 12.01817286, 3.8537762, 5.6709554, 7.32330843, 5.12760355, 15.59984166, 16.02646317, 10.17028229, 8.15392771, 17.89921969, 1.45500612, 6.04264206, 13.28824345, 16.18839369, 11.69951924, 6.08602176] ,[12.64895326, 5.95639905, 1.75205743, 6.15052067, 8.17549729, 9.85930991, 7.93540704, 11.55480272, 3.69386905, 4.58778427, 17.44387647, 11.42922451, 12.92102814, 0.0, 4.27095536, 0.94275412, 12.86737946, 10.79908569, 5.98328464, 14.06290555, 12.00703063, 11.18734576, 2.85685962, 9.05534815, 10.94566708, 12.13076671, 15.99908456, 0.92611969, 9.94017355, 1.66376675, 13.59119982] ,[11.83495959, 2.40239869, 2.88335684, 2.34510263, 4.66479455, 6.81434145, 4.63742738, 9.05636221, 5.70986429, 5.46866104, 14.25705759, 8.8175326, 11.20494726, 4.27095536, 0.0, 4.10243882, 9.9869739, 7.50100139, 6.22658749, 10.89604224, 8.12961185, 7.56920751, 4.30290909, 5.18641268, 8.12049236, 10.05833251, 12.81750885, 3.80408544, 6.70466442, 2.72495174, 10.1660557] ,[11.70879828, 6.07824264, 2.16997718, 5.70239291, 7.58322465, 9.13965091, 8.17827032, 12.03251743, 2.81169309, 3.6463882, 16.69623841, 10.62470016, 12.01817286, 0.94275412, 4.10243882, 0.0, 12.09053409, 10.11629031, 5.04063697, 13.32719632, 12.10124288, 11.37477017, 1.91520346, 8.54783656, 11.34642484, 11.26028144, 15.25622608, 1.67690854, 10.22931137, 1.38448682, 12.91187247] ,[6.4373446, 11.11678063, 12.37288352, 7.69582198, 5.41598382, 3.17264835, 12.37857421, 16.55079062, 10.9383549, 9.54749941, 4.65581362, 1.58412805, 3.8537762, 12.86737946, 9.9869739, 12.09053409, 0.0, 2.71236511, 8.45717468, 1.54718684, 12.49655477, 13.14306389, 10.52319648, 5.55523918, 15.32051453, 2.62909282, 3.26225878, 12.96975405, 13.61170773, 11.35459032, 2.24101784] ,[7.78414476, 8.44551052, 10.07524752, 5.1626248, 2.83912446, 1.12205638, 9.66627221, 13.86878021, 9.48981659, 8.1910259, 6.76328816, 2.34615336, 5.6709554, 10.79908569, 7.50100139, 10.11629031, 2.71236511, 0.0, 7.4137736, 3.39538331, 9.99615576, 10.54206128, 8.78763952, 2.84316438, 12.64364243, 4.21652771, 5.31945246, 10.77035217, 10.93028738, 9.19749848, 2.79562155] ,[6.68451163, 8.60010133, 6.58401409, 5.78125939, 6.06284458, 6.29177839, 10.80194322, 15.28080044, 2.75163283, 1.43237716, 13.06730703, 6.87329623, 7.32330843, 5.98328464, 6.22658749, 5.04063697, 8.45717468, 7.4137736, 0.0, 9.93141536, 13.66615417, 13.428634, 3.12823272, 7.36493564, 14.31770932, 6.88185392, 11.71322338, 6.58013675, 12.8215552, 5.22461902, 9.95537903] ,[7.79300168, 11.75839269, 13.4353958, 8.55787492, 6.23423946, 4.2081104, 12.78714438, 16.74021478, 12.34462784, 10.96807515, 3.39012822, 3.09159006, 5.12760355, 14.06290555, 10.89604224, 13.32719632, 1.54718684, 3.39538331, 9.93141536, 0.0, 12.37570065, 13.16726412, 11.84428103, 6.05968734, 15.50032575, 4.07283787, 1.94123626, 14.09525309, 13.8187458, 12.49965066, 1.29839868] ,[17.34996024, 6.05118257, 10.28382056, 7.88908566, 8.50188127, 10.28569975, 4.30212228, 5.10872591, 13.79041051, 13.31152777, 14.82726971, 12.31382888, 15.59984166, 12.00703063, 8.12961185, 12.10124288, 12.49655477, 9.99615576, 13.66615417, 12.37570065, 0.0, 1.39939368, 12.40208116, 7.50154126, 4.01420874, 14.15266741, 13.72428755, 11.25216342, 3.05096951, 10.72886707, 11.11246763] ,[17.60395722, 5.30404834, 9.44039185, 7.73356999, 8.71843188, 10.69209409, 3.27891144, 3.83728518, 13.27898131, 12.93395999, 15.79334966, 12.80286736, 16.02646317, 11.18734576, 7.56920751, 11.37477017, 13.14306389, 10.54206128, 13.428634, 13.16726412, 1.39939368, 0.0, 11.87174275, 7.87263384, 2.67596747, 14.59820309, 14.62670887, 10.38281294, 1.68776206, 10.03049275, 11.93787468] ,[9.79397585, 6.63675716, 3.62189896, 5.10321362, 6.53196685, 7.73901619, 8.86209542, 13.06943212, 1.41196102, 1.75956695, 15.16990059, 9.00424136, 10.17028229, 2.85685962, 4.30290909, 1.91520346, 10.52319648, 8.78763952, 3.12823272, 11.84428103, 12.40208116, 11.87174275, 0.0, 7.66463837, 12.2496542, 9.48467094, 13.74649513, 3.47565239, 10.94402761, 2.31264759, 11.55864464] ,[9.85718622, 5.70534735, 7.99620437, 2.92762016, 1.30475623, 2.81956828, 6.82343432, 11.06218161, 8.70828004, 7.65194891, 9.2924934, 4.95939238, 8.15392771, 9.05534815, 5.18641268, 8.54783656, 5.55523918, 2.84316438, 7.36493564, 6.05968734, 7.50154126, 7.87263384, 7.66463837, 0.0, 9.84604144, 6.72738767, 7.89198533, 8.82474647, 8.12925592, 7.39221693, 5.12787027] ,[19.21619553, 5.72297753, 9.22408933, 8.99537615, 10.46645486, 12.62381157, 3.52381635, 1.24056048, 13.65033066, 13.57881789, 18.2779851, 14.7979951, 17.89921969, 10.94566708, 8.12049236, 11.34642484, 15.32051453, 12.64364243, 14.31770932, 15.50032575, 4.01420874, 2.67596747, 12.2496542, 9.84604144, 0.0, 16.51298451, 17.05842325, 10.05243704, 1.71703268, 10.12636373, 14.3150126] ,[3.82862739, 11.6685413, 11.99509716, 8.01345658, 6.07234356, 3.90893715, 13.28684489, 17.70810142, 9.5787741, 8.17796847, 6.51126118, 1.87245486, 1.45500612, 12.13076671, 10.05833251, 11.26028144, 2.62909282, 4.21652771, 6.88185392, 4.07283787, 14.15266741, 14.59820309, 9.48467094, 6.72738767, 16.51298451, 0.0, 5.34972425, 12.41639362, 14.79886643, 10.80266929, 4.85890174] ,[8.7164367, 13.59625527, 15.37513181, 10.47528224, 8.15276462, 6.14837579, 14.50755284, 18.29819095, 14.19733006, 12.80860634, 1.44894598, 4.84219825, 6.04264206, 15.99908456, 12.81750885, 15.25622608, 3.26225878, 5.31945246, 11.71322338, 1.94123626, 13.72428755, 14.62670887, 13.74649513, 7.89198533, 17.05842325, 5.34972425, 0.0, 16.03638109, 15.4070327, 14.4399125, 2.79937419] ,[13.18288407, 5.22122152, 0.98681961, 5.89713902, 8.05259831, 9.88932377, 7.11414623, 10.64149624, 4.48787087, 5.23515614, 17.48502524, 11.58558897, 13.28824345, 0.92611969, 3.80408544, 1.67690854, 12.96975405, 10.77035217, 6.58013675, 14.09525309, 11.25216342, 10.38281294, 3.47565239, 8.82474647, 10.05243704, 12.41639362, 16.03638109, 0.0, 9.09019728, 1.624029, 13.54359517] ,[17.53956817, 4.31369829, 8.18963269, 7.36247894, 8.75921413, 10.90755626, 2.08216715, 2.93918273, 12.35512155, 12.16922105, 16.64965172, 13.08101108, 16.18839369, 9.94017355, 6.70466442, 10.22931137, 13.61170773, 10.93028738, 12.8215552, 13.8187458, 3.05096951, 1.68776206, 10.94402761, 8.12925592, 1.71703268, 14.79886643, 15.4070327, 9.09019728, 0.0, 8.93683876, 12.64820729] ,[11.70369435, 4.7267962, 1.3594828, 4.48740453, 6.53127564, 8.29159556, 6.86896002, 10.87722827, 3.61488034, 3.99564885, 15.88775764, 9.96271172, 11.69951924, 1.66376675, 2.72495174, 1.38448682, 11.35459032, 9.19749848, 5.22461902, 12.49965066, 10.72886707, 10.03049275, 2.31264759, 7.39221693, 10.12636373, 10.80266929, 14.4399125, 1.624029, 8.93683876, 0.0, 11.98209726] ,[8.67701736, 10.82099606, 12.82301273, 7.82270527, 5.52542753, 3.82122721, 11.70959286, 15.55557199, 12.19699495, 10.86110176, 4.16856321, 3.44220958, 6.08602176, 13.59119982, 10.1660557, 12.91187247, 2.24101784, 2.79562155, 9.95537903, 1.29839868, 11.11246763, 11.93787468, 11.55864464, 5.12787027, 14.3150126, 4.85890174, 2.79937419, 13.54359517, 12.64820729, 11.98209726, 0.0]
    # AdjList_2
    mov = [0.0, 9.1822057, 5.37870005, 6.47364319, 9.53819148, 5.92046594, 4.88462155, 6.52737148, 11.39766237, 5.74520131, 4.71163119, 14.44962547, 6.19092425, 5.54268602, 10.1701806, 9.73094642, 4.00147642, 5.54034297, 8.06971741, 2.6622319, 3.0316653, 6.17529954, 7.53095923, 6.47593835, 9.57021193, 7.02924044, 13.15615481, 1.4155841, 1.91474863, 1.03061884, 4.08976547] ,[9.1822057, 0.0, 10.85766452, 2.70867108, 8.32843247, 4.36889397, 4.97211509, 14.31688529, 2.23099713, 12.68860205, 9.73125231, 5.38429567, 9.8239692, 8.81836576, 1.45213432, 1.35025323, 10.9632914, 9.14376248, 7.40499404, 11.70086796, 9.46478606, 6.63378412, 15.68492473, 8.78894116, 2.38736413, 10.66970362, 4.03169141, 10.26538179, 11.09049126, 9.94767034, 7.55924726] ,[5.37870005, 10.85766452, 0.0, 8.59621157, 6.41855433, 6.49289212, 8.62218046, 4.0191827, 12.70648999, 2.13037585, 10.08017726, 15.16495173, 11.5204844, 2.29257495, 12.22678194, 10.71843068, 9.26235094, 10.82213547, 5.46892689, 6.72432892, 8.40983227, 4.61546109, 5.48421765, 2.96875411, 9.97855428, 12.39031392, 14.13921229, 6.38376656, 5.45670572, 4.68028726, 9.08432841] ,[6.47364319, 2.70867108, 8.59621157, 0.0, 7.60889501, 2.48382385, 2.64698588, 11.80432445, 4.92782759, 10.27121289, 7.48116129, 8.02732687, 7.87774355, 6.7928557, 3.8090798, 3.41058872, 8.49027577, 7.11903194, 6.33097869, 9.0066879, 6.96113036, 4.95571532, 13.13173667, 7.00400059, 3.66118737, 8.80823595, 6.7044447, 7.57559678, 8.38193656, 7.24678416, 5.2745954] ,[9.53819148, 8.32843247, 6.41855433, 7.60889501, 0.0, 5.31723634, 9.46765217, 10.36098563, 9.21916514, 8.48343844, 13.29163879, 10.61214746, 14.29793045, 4.46365382, 9.73587393, 7.42376573, 13.34166634, 13.4912253, 1.51595183, 11.89250316, 11.98541149, 3.475449, 11.79682063, 3.53319268, 6.27700821, 15.27807571, 10.01109445, 10.93522684, 10.67789962, 9.46914919, 11.42568282] ,[5.92046594, 4.36889397, 6.49289212, 2.48382385, 5.31723634, 0.0, 4.18150555, 10.06997013, 6.25997574, 8.37270422, 8.46079745, 8.96569997, 9.2419498, 4.47279551, 5.76269941, 4.31025857, 8.97269349, 8.43748865, 3.935713, 8.58222177, 7.46659869, 2.4759381, 11.47899848, 4.56559744, 3.80887614, 10.21879529, 7.81102629, 7.26609121, 7.66505454, 6.39573712, 6.39220047] ,[4.88462155, 4.97211509, 8.62218046, 2.64698588, 9.46765217, 4.18150555, 0.0, 11.08359687, 7.15993207, 9.86023014, 4.84624339, 10.3418945, 5.25971621, 7.3904523, 5.60496406, 5.95045856, 5.99159466, 4.48676534, 8.02641835, 7.06989926, 4.50253828, 6.25752578, 12.25149515, 7.91769317, 6.30612915, 6.20827338, 8.97333669, 5.64333453, 6.7435628, 5.84936448, 2.64007794] ,[6.52737148, 14.31688529, 4.0191827, 11.80432445, 10.36098563, 10.06997013, 11.08359687, 0.0, 16.32379008, 1.89193962, 10.79445649, 18.96839592, 12.35825555, 6.28308351, 15.58997329, 14.37897522, 9.35900749, 11.83428359, 9.48525717, 6.18596159, 9.14197915, 8.53589682, 1.46760229, 6.98321615, 13.76010436, 13.03968822, 17.8676216, 6.73234295, 5.32458869, 5.50183654, 10.60808239] ,[11.39766237, 2.23099713, 12.70648999, 4.92782759, 9.21916514, 6.25997574, 7.15993207, 16.32379008, 0.0, 14.63053723, 11.83971748, 3.1857805, 11.80420214, 10.5526888, 1.90883365, 1.99939637, 13.14647093, 11.16805764, 8.60254208, 13.92966956, 11.66201398, 8.25197814, 17.72083036, 10.37205992, 2.98324778, 12.59823958, 1.81341764, 12.4950742, 13.30138195, 12.14087177, 9.70084124] ,[5.74520131, 12.68860205, 2.13037585, 10.27121289, 8.48343844, 8.37270422, 9.86023014, 1.89193962, 14.63053723, 0.0, 10.3482718, 17.19278267, 11.8795482, 4.40862822, 14.00799259, 12.66252255, 9.17423412, 11.26750584, 7.59352594, 6.20324606, 8.65129624, 6.69426301, 3.35445164, 5.09141403, 11.98877583, 12.65749028, 16.12602131, 6.33573863, 5.07890873, 4.80036522, 9.78847366] ,[4.71163119, 9.73125231, 10.08017726, 7.48116129, 13.29163879, 8.46079745, 4.84624339, 10.79445649, 11.83971748, 10.3482718, 0.0, 15.01785958, 1.56555349, 9.92255276, 10.10868123, 10.78456743, 1.90813972, 1.22483434, 11.7761492, 4.83686291, 1.69864353, 9.82193319, 11.49057053, 10.76037116, 11.14156574, 2.31983302, 13.63749401, 4.06857557, 5.47102029, 5.56004406, 2.20717077] ,[14.44962547, 5.38429567, 15.16495173, 8.02732687, 10.61214746, 8.96569997, 10.3418945, 18.96839592, 3.1857805, 17.19278267, 15.01785958, 0.0, 14.93449838, 12.91183809, 4.95917399, 4.73557379, 16.33081353, 14.31926175, 10.41035252, 17.0258604, 14.84236711, 10.56396584, 20.40198242, 12.54892626, 5.20830942, 15.6997856, 1.38071457, 15.5994272, 16.33479482, 15.13665846, 12.88573596] ,[6.19092425, 9.8239692, 11.5204844, 7.87774355, 14.29793045, 9.2419498, 5.25971621, 12.35825555, 11.80420214, 11.8795482, 1.56555349, 14.93449838, 0.0, 11.1993975, 9.97585752, 10.99250031, 3.3158576, 0.80796721, 12.79227355, 6.38321177, 3.2355551, 10.86507049, 13.05317708, 11.99066116, 11.50183455, 0.98109287, 13.55902137, 5.6290112, 7.03577417, 7.08126105, 2.87438739] ,[5.54268602, 8.81836576, 2.29257495, 6.7928557, 4.46365382, 4.47279551, 7.3904523, 6.28308351, 10.5526888, 4.40862822, 9.92255276, 12.91183809, 11.1993975, 0.0, 10.23032355, 8.55440014, 9.53229199, 10.43256283, 3.29483546, 7.60933351, 8.38099472, 2.34792997, 7.75068104, 0.97984375, 7.75114213, 12.14012532, 11.92123466, 6.84599979, 6.35184578, 5.23264882, 8.48394921] ,[10.1701806, 1.45213432, 12.22678194, 3.8090798, 9.73587393, 5.76269941, 5.60496406, 15.58997329, 1.90883365, 14.00799259, 10.10868123, 4.95917399, 9.97585752, 10.23032355, 0.0, 2.42969943, 11.51960264, 9.36874926, 8.84894732, 12.58750803, 10.07681834, 8.07396239, 16.93199251, 10.22834858, 3.60352933, 10.74347513, 3.58707588, 11.14754206, 12.08461203, 11.00233599, 8.02565824] ,[9.73094642, 1.35025323, 10.71843068, 3.41058872, 7.42376573, 4.31025857, 5.95045856, 14.37897522, 1.99939637, 12.66252255, 10.78456743, 4.73557379, 10.99250031, 8.55440014, 2.42969943, 0.0, 11.8888578, 10.2835447, 6.68342381, 12.33139264, 10.36469279, 6.2531223, 15.78911974, 8.37743876, 1.18127017, 11.86662168, 3.51236977, 10.91619771, 11.60673442, 10.40212215, 8.58515708] ,[4.00147642, 10.9632914, 9.26235094, 8.49027577, 13.34166634, 8.97269349, 5.99159466, 9.35900749, 13.14647093, 9.17423412, 1.90813972, 16.33081353, 3.3158576, 9.53229199, 11.51960264, 11.8888578, 0.0, 3.130995, 11.83803171, 3.20033242, 1.53764459, 9.88691963, 9.90592446, 10.44926213, 12.08461827, 3.77509409, 14.95919007, 2.89019543, 4.10885803, 4.58306233, 3.51306396] ,[5.54034297, 9.14376248, 10.82213547, 7.11903194, 13.4912253, 8.43748865, 4.48676534, 11.83428359, 11.16805764, 11.26750584, 1.22483434, 14.31926175, 0.80796721, 10.43256283, 9.36874926, 10.2835447, 3.130995, 0.0, 11.98626089, 6.00040854, 2.69857762, 10.06127539, 12.59480662, 11.21155597, 10.75811181, 1.78904269, 12.94024844, 5.10914533, 6.53478795, 6.46932895, 2.06664059] ,[8.06971741, 7.40499404, 5.46892689, 6.33097869, 1.51595183, 3.935713, 8.02641835, 9.48525717, 8.60254208, 7.59352594, 11.7761492, 10.41035252, 12.79227355, 3.29483546, 8.84894732, 6.68342381, 11.83803171, 11.98626089, 0.0, 10.47967448, 10.47288191, 1.96041786, 10.94679423, 2.502159, 5.62010379, 13.7718087, 9.63418793, 9.47502145, 9.28906148, 8.05504399, 9.92207298] ,[2.6622319, 11.70086796, 6.72432892, 9.0066879, 11.89250316, 8.58222177, 7.06989926, 6.18596159, 13.92966956, 6.20324606, 4.83686291, 17.0258604, 6.38321177, 7.60933351, 12.58750803, 12.33139264, 3.20033242, 6.00040854, 10.47967448, 0.0, 3.3843097, 8.66419348, 6.70856393, 8.58664694, 12.2169301, 6.95350444, 15.71100705, 1.44024793, 1.28000536, 2.42503482, 5.40414745] ,[3.0316653, 9.46478606, 8.40983227, 6.96113036, 11.98541149, 7.46659869, 4.50253828, 9.14197915, 11.66201398, 8.65129624, 1.69864353, 14.84236711, 3.2355551, 8.38099472, 10.07681834, 10.36469279, 1.53764459, 2.69857762, 10.47288191, 3.3843097, 0.0, 8.51324761, 9.89939246, 9.25620956, 10.54758964, 4.0163829, 13.4754256, 2.41169773, 3.83621367, 3.86142483, 2.19658771] ,[6.17529954, 6.63378412, 4.61546109, 4.95571532, 3.475449, 2.4759381, 6.25752578, 8.53589682, 8.25197814, 6.69426301, 9.82193319, 10.56396584, 10.86507049, 2.34792997, 8.07396239, 6.2531223, 9.88691963, 10.06127539, 1.96041786, 8.66419348, 8.51324761, 0.0, 9.99806396, 2.15517825, 5.41008357, 11.84260452, 9.58070977, 7.5888626, 7.5199737, 6.25390261, 8.00256223] ,[7.53095923, 15.68492473, 5.48421765, 13.13173667, 11.79682063, 11.47899848, 12.25149515, 1.46760229, 17.72083036, 3.35445164, 11.49057053, 20.40198242, 13.05317708, 7.75068104, 16.93199251, 15.78911974, 9.90592446, 12.59480662, 10.94679423, 6.70856393, 9.89939246, 9.99806396, 0.0, 8.4446718, 15.19483111, 13.66012198, 19.28634271, 7.50259225, 6.07322933, 6.5026104, 11.545927] ,[6.47593835, 8.78894116, 2.96875411, 7.00400059, 3.53319268, 4.56559744, 7.91769317, 6.98321615, 10.37205992, 5.09141403, 10.76037116, 12.54892626, 11.99066116, 0.97984375, 10.22834858, 8.37743876, 10.44926213, 11.21155597, 2.502159, 8.58664694, 9.25620956, 2.15517825, 8.4446718, 0.0, 7.47715745, 12.94349532, 11.63369509, 7.80004064, 7.33094153, 6.20287638, 9.22329968] ,[9.57021193, 2.38736413, 9.97855428, 3.66118737, 6.27700821, 3.80887614, 6.30612915, 13.76010436, 2.98324778, 11.98877583, 11.14156574, 5.20830942, 11.50183455, 7.75114213, 3.60352933, 1.18127017, 12.08461827, 10.75811181, 5.62010379, 12.2169301, 10.54758964, 5.41008357, 15.19483111, 7.47715745, 0.0, 12.41263908, 4.17068181, 10.836299, 11.39202245, 10.1456389, 8.93546999] ,[7.02924044, 10.66970362, 12.39031392, 8.80823595, 15.27807571, 10.21879529, 6.20827338, 13.03968822, 12.59823958, 12.65749028, 2.31983302, 15.6997856, 0.98109287, 12.14012532, 10.74347513, 11.86662168, 3.77509409, 1.78904269, 13.7718087, 6.95350444, 4.0163829, 11.84260452, 13.66012198, 12.94349532, 12.41263908, 0.0, 14.33052296, 6.34748161, 7.72041625, 7.8768429, 3.85543646] ,[13.15615481, 4.03169141, 14.13921229, 6.7044447, 10.01109445, 7.81102629, 8.97333669, 17.8676216, 1.81341764, 16.12602131, 13.63749401, 1.38071457, 13.55902137, 11.92123466, 3.58707588, 3.51236977, 14.95919007, 12.94024844, 9.63418793, 15.71100705, 13.4754256, 9.58070977, 19.28634271, 11.63369509, 4.17068181, 14.33052296, 0.0, 14.27980435, 15.05156559, 13.87084461, 11.50803563] ,[1.4155841, 10.26538179, 6.38376656, 7.57559678, 10.93522684, 7.26609121, 5.64333453, 6.73234295, 12.4950742, 6.33573863, 4.06857557, 15.5994272, 5.6290112, 6.84599979, 11.14754206, 10.91619771, 2.89019543, 5.10914533, 9.47502145, 1.44024793, 2.41169773, 7.5888626, 7.50259225, 7.80004064, 10.836299, 6.34748161, 14.27980435, 0.0, 1.42972012, 1.70668492, 4.15974727] ,[1.91474863, 11.09049126, 5.45670572, 8.38193656, 10.67789962, 7.66505454, 6.7435628, 5.32458869, 13.30138195, 5.07890873, 5.47102029, 16.33479482, 7.03577417, 6.35184578, 12.08461203, 11.60673442, 4.10885803, 6.53478795, 9.28906148, 1.28000536, 3.83621367, 7.5199737, 6.07322933, 7.33094153, 11.39202245, 7.72041625, 15.05156559, 1.42972012, 0.0, 1.28838946, 5.53672994] ,[1.03061884, 9.94767034, 4.68028726, 7.24678416, 9.46914919, 6.39573712, 5.84936448, 5.50183654, 12.14087177, 4.80036522, 5.56004406, 15.13665846, 7.08126105, 5.23264882, 11.00233599, 10.40212215, 4.58306233, 6.46932895, 8.05504399, 2.42503482, 3.86142483, 6.25390261, 6.5026104, 6.20287638, 10.1456389, 7.8768429, 13.87084461, 1.70668492, 1.28838946, 0.0, 5.10764123] ,[4.08976547, 7.55924726, 9.08432841, 5.2745954, 11.42568282, 6.39220047, 2.64007794, 10.60808239, 9.70084124, 9.78847366, 2.20717077, 12.88573596, 2.87438739, 8.48394921, 8.02565824, 8.58515708, 3.51306396, 2.06664059, 9.92207298, 5.40414745, 2.19658771, 8.00256223, 11.545927, 9.22329968, 8.93546999, 3.85543646, 11.50803563, 4.15974727, 5.53672994, 5.10764123, 0.0]
    # AdjList_3
    #mov = [[0.0, 16.3373068, 12.59291712, 4.63741258, 1.97803401, 8.52592732, 6.51121033, 0.94572927, 14.79800582, 15.83409495, 1.01541926, 13.54466405, 9.03049571, 10.4, 17.68, 15.03, 16.68, 16.21, 18.12, 2.79, 18.76] ,[16.3373068, 0.0, 4.54843225, 12.36515352, 14.55772824, 7.98274869, 10.25717933, 17.10224018, 2.26981913, 3.25244701, 15.47551663, 6.05800402, 8.21332868, 5.97029581, 1.94018838, 7.75940753, 9.31749967, 4.73548484, 10.58588664, 14.19190572, 3.47467575] ,[12.59291712, 4.54843225, 0.0, 9.31083611, 10.68434463, 5.23702399, 7.38323237, 13.44484348, 2.41935334, 3.25327204, 11.84907791, 2.57085697, 6.28390112, 2.65524993, 6.41846523, 4.93282399, 6.96472041, 3.95918532, 8.55770996, 10.87826423, 7.90053415] ,[4.63741258, 12.36515352, 9.31083611, 0.0, 3.72136656, 4.38668941, 2.11204417, 5.12892403, 11.18137531, 12.51254302, 3.62582112, 10.95104198, 4.50417707, 6.75597478, 13.49782143, 12.93103026, 14.86189778, 13.22993455, 16.45192709, 1.94654952, 14.43516664] ,[1.97803401, 14.55772824, 10.68434463, 3.72136656, 0.0, 6.93545464, 5.17573595, 2.90817787, 12.93341212, 13.91265403, 1.63001742, 11.57096567, 7.67800117, 8.59190439, 15.98891396, 13.06173108, 14.72787817, 14.25880998, 16.18383983, 2.52203503, 17.12489712] ,[8.52592732, 7.98274869, 5.23702399, 4.38668941, 6.93545464, 0.0, 2.29182564, 9.21920738, 6.85475821, 8.29651474, 7.60144426, 7.33370519, 1.5774714, 2.58177444, 9.19123479, 9.58974767, 11.63998379, 9.19165469, 13.27302084, 6.2180581, 10.23460197] ,[6.51121033, 10.25717933, 7.38323237, 2.11204417, 5.17573595, 2.29182564, 0.0, 7.11589522, 9.13268993, 10.52751953, 5.5343602, 9.24118711, 2.5258809, 4.76270493, 11.39113903, 11.3643853, 13.36479887, 11.33683809, 14.98438175, 4.00607307, 12.35188631] ,[0.94572927, 17.10224018, 13.44484348, 5.12892403, 2.90817787, 9.21920738, 7.11589522, 0.0, 15.61759823, 16.69255143, 1.62674416, 14.45541487, 9.60745105, 11.19663129, 18.40701744, 15.96984174, 17.62496808, 17.10740229, 19.06665787, 3.19439148, 19.43912997] ,[14.79800582, 2.26981913, 2.41935334, 11.18137531, 12.93341212, 6.85475821, 9.13268993, 15.61759823, 0.0, 1.75759664, 14.00180368, 3.82625933, 7.49723056, 4.45530971, 4.20759192, 5.72253428, 7.47297927, 3.1878841, 8.88708292, 12.88654043, 5.73769882] ,[15.83409495, 3.25244701, 3.25327204, 12.51254302, 13.91265403, 8.29651474, 10.52751953, 16.69255143, 1.75759664, 0.0, 15.10029259, 3.30845534, 9.09422052, 5.76702633, 5.05629969, 4.57399843, 6.06509665, 1.51883639, 7.36066297, 14.12088785, 6.56260814] ,[1.01541926, 15.47551663, 11.84907791, 3.62582112, 1.63001742, 7.60144426, 5.5343602, 1.62674416, 14.00180368, 15.10029259, 0.0, 12.94127508, 8.04529259, 9.57430267, 16.78544578, 14.54105786, 16.25783837, 15.55176795, 17.73716703, 1.78834779, 17.8288574] ,[13.54466405, 6.05800402, 2.57085697, 10.95104198, 11.57096567, 7.33370519, 9.24118711, 14.45541487, 3.82625933, 3.30845534, 12.94127508, 0.0, 8.58980729, 4.86573565, 7.99379521, 2.37450623, 4.42911525, 2.87482559, 6.04499953, 12.26887452, 9.53036185] ,[9.03049571, 8.21332868, 6.28390112, 4.50417707, 7.67800117, 1.5774714, 2.5258809, 9.60745105, 7.49723056, 9.09422052, 8.04529259, 8.58980729, 0.0, 3.72432256, 9.12509814, 10.9090564, 12.97251228, 10.15495101, 14.60370641, 6.44838888, 9.97119374] ,[10.40151068, 5.97029581, 2.65524993, 6.75597478, 8.59190439, 2.58177444, 4.76270493, 11.19663129, 4.45530971, 5.76702633, 9.57430267, 4.86573565, 3.72432256, 0.0, 7.5015882, 7.19536555, 9.2607303, 6.61118705, 10.88943006, 8.43170198, 8.77548441] ,[17.68984872, 1.94018838, 6.41846523, 13.49782143, 15.98891396, 9.19123479, 11.39113903, 18.40701744, 4.20759192, 5.05629969, 16.78544578, 7.99379521, 9.12509814, 7.5015882, 0.0, 9.6207048, 11.07533017, 6.46311062, 12.25236496, 15.3884362, 1.53667181] ,[15.03443016, 7.75940753, 4.93282399, 12.93103026, 13.06173108, 9.58974767, 11.3643853, 15.96984174, 5.72253428, 4.57399843, 14.54105786, 2.37450623, 10.9090564, 7.19536555, 9.6207048, 0.0, 2.06563516, 3.35163053, 3.69474228, 14.07701588, 11.13479968] ,[16.68166788, 9.31749967, 6.96472041, 14.86189778, 14.72787817, 11.63998379, 13.36479887, 17.62496808, 7.47297927, 6.06509665, 16.25783837, 4.42911525, 12.97251228, 9.2607303, 11.07533017, 2.06563516, 0.0, 4.61222109, 1.6332179, 15.9141725, 12.54030497] ,[16.2170464, 4.73548484, 3.95918532, 13.22993455, 14.25880998, 9.19165469, 11.33683809, 17.10740229, 3.1878841, 1.51883639, 15.55176795, 2.87482559, 10.15495101, 6.61118705, 6.46311062, 3.35163053, 4.61222109, 0.0, 5.8529282, 14.72013023, 7.93355686] ,[18.12131512, 10.58588664, 8.55770996, 16.45192709, 16.18383983, 13.27302084, 14.98438175, 19.06665787, 8.88708292, 7.36066297, 17.73716703, 6.04499953, 14.60370641, 10.88943006, 12.25236496, 3.69474228, 1.6332179, 5.8529282, 0.0, 17.45774002, 13.66747856] ,[2.79602033, 14.19190572, 10.87826423, 1.94654952, 2.52203503, 6.2180581, 4.00607307, 3.19439148, 12.88654043, 14.12088785, 1.78834779, 12.26887452, 6.44838888, 8.43170198, 15.3884362, 14.07701588, 15.9141725, 14.72013023, 17.45774002, 0.0, 16.35707464] ,[18.76046551, 3.47467575, 7.90053415, 14.43516664, 17.12489712, 10.23460197, 12.35188631, 19.43912997, 5.73769882, 6.56260814, 17.8288574, 9.53036185, 9.97119374, 8.77548441, 1.53667181, 11.13479968, 12.54030497, 7.93355686, 13.66747856, 16.35707464, 0.0] ]
    #T = 32 # rondas de defensa (debe ser mayor o igual al número máximo de defendidos)
    T2 = 12 # rondas de quema locales
    Tb = 20 # rondas de quema globales
    best_defense = []
    best_defense_time = []
    total_runtime = 0
    best_solution = []
    size_best_solution = float("inf")
    for i in range(len(dataset)):
        file = folder_dataset + dataset[i][0]
        g = nx.read_adjlist(file, nodetype=int)
        G = g.copy()
        print("--------------------------------------------------------------")
        instance = dataset[i][0]
        print("instance: " + instance)
        main(len(G.nodes()), len(G.edges()), folder_dataset + dataset[i][0], dataset[i][1])
        upper = 20
        lower = 1
        while lower <= upper:
            mid = math.floor((upper + lower) /2)
            T = mid # Aquí ya sabemos cuantas rondas de defensa hay
            print("MID: " + str(mid))
            run(T, T2, mid, Tb)
            total_runtime = total_runtime + runtime
            print("is feasible: " + str(feasible))
            
            if not model_feasible:
                upper = mid - 1
            else:
                if feasible:
                    upper = mid - 1
                else:
                    lower = mid + 1
        print("defense: " + str(best_solution[0]))
        print("defense time: " + str(best_solution[1]))
        print("time: " + str(total_runtime))
