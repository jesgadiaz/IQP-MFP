from gurobipy import *
import math
import numpy as np
import time
import networkx as nx
import random
import copy

def createGraph(input_file):
    global G, n, m, a0, start_time, apsp_time, n_nodes, n_edges, conn_comp, av_degree, density, card_comp, delta_m
    
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

def run(T_input, T2, Td):
    global G, total_runtime, runtime, n, m, a0, feasible, best_sequence, I_input, delta_m, feasible, best_defense, best_defense_time
    try:
        m = Model("mip1")
        m.Params.outputFlag = 1  # 0 - Off  //  1 - On
        m.setParam("MIPGap", 0.0);
        m.setParam("Presolve", 2); # -1 - Automatic // 0 - Off // 1 - Conservative // 2 - Aggresive 
        #m.setParam("PreQLinearize", -1); # -1 - Automatic // 0 - Off // 1 - Strong LP relaxation // 2 - Compact relaxation
        #m.params.BestObjStop = k
        
        T = T_input
        M = n + 1
        
        
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
            
        mov = [[0.0, 2.69453808, 3.64161647, 12.34325843, 2.86009481, 6.1374802, 5.32246651, 11.29659716, 3.31525444, 0.78519918, 12.41331877, 12.71590282, 10.85152542, 5.30531019, 14.62988656, 15.66729075, 5.1738933, 3.58142053, 8.10023862, 7.72959113, 2.08169749, 4.59862711, 6.84565681, 6.91266051, 16.3123257, 4.8405079], [2.69453808, 0.0, 6.33286266, 9.64875876, 4.29425366, 3.62836634, 6.02264541, 8.60592167, 3.06480263, 2.25476522, 9.73701553, 10.12917852, 8.17368452, 7.60734617, 12.13698276, 13.27478887, 2.47957666, 5.94352464, 5.43828894, 5.04438059, 4.76417279, 2.78090837, 7.02480049, 7.70310563, 13.81340365, 6.77107657], [3.64161647, 6.33286266, 0.0, 15.9780014, 4.00804448, 9.64134273, 6.17605611, 14.93774382, 5.90190981, 4.16710899, 16.01513254, 16.21524974, 14.45748061, 3.5593164, 18.01879751, 18.94081739, 8.81211911, 2.56523183, 11.73700783, 11.34952072, 1.57468355, 7.8158263, 7.95886858, 7.21101714, 19.69790598, 4.25692634], [12.34325843, 9.64875876, 15.9780014, 0.0, 13.18906052, 6.65734739, 13.53490248, 1.21019963, 11.06493431, 11.83677242, 1.21847207, 3.11832915, 1.83175951, 16.84611761, 5.18439858, 6.82614445, 7.17018836, 15.33762043, 4.44343542, 4.65384074, 14.40544845, 8.8601632, 13.28726193, 14.79876271, 6.2995567, 15.67399502], [2.86009481, 4.29425366, 4.00804448, 13.18906052, 0.0, 6.54302706, 2.55203149, 12.30111758, 2.2016778, 2.46256654, 12.99438637, 12.87478936, 11.49295696, 7.07799438, 14.46175496, 15.24009862, 6.40054883, 5.49462281, 9.40801035, 8.57909014, 2.92421481, 4.40130551, 4.23596374, 4.06829594, 16.11947556, 7.11223997], [6.1374802, 3.62836634, 9.64134273, 6.65734739, 6.54302706, 0.0, 7.02023245, 5.83650559, 4.40798494, 5.49154145, 6.46790425, 6.58377272, 4.9515114, 11.23322217, 8.51841325, 9.64893718, 1.89886861, 9.55770251, 3.47662582, 2.17386813, 8.07913676, 2.21692624, 7.13140097, 8.44665983, 10.19802227, 10.37950982], [5.32246651, 6.02264541, 6.17605611, 13.53490248, 2.55203149, 7.02023245, 0.0, 12.82257679, 3.01623754, 4.77488592, 13.11344175, 12.62301977, 11.7338606, 9.50479511, 13.85853526, 14.38249205, 7.54307763, 7.98929774, 10.34049008, 9.1927582, 5.38282355, 4.84396878, 1.78281258, 1.68452975, 15.44610294, 9.64217133], [11.29659716, 8.60592167, 14.93774382, 1.21019963, 12.30111758, 5.83650559, 12.82257679, 0.0, 10.2247338, 10.82382343, 2.09441163, 3.84432566, 1.64695079, 15.6963289, 6.09367811, 7.72561094, 6.12954017, 14.21040127, 3.28797438, 3.72233099, 13.37009259, 8.05330675, 12.71380328, 14.16409044, 7.34540194, 14.50372216], [3.31525444, 3.06480263, 5.90190981, 11.06493431, 2.2016778, 4.40798494, 3.01623754, 10.2247338, 0.0, 2.55091774, 10.8229774, 10.67314007, 9.33868012, 8.46980964, 12.27957352, 13.1003455, 4.58428986, 6.76531945, 7.4908784, 6.51112833, 4.51501243, 2.21931093, 3.98174339, 4.67500764, 13.94327444, 8.14568881], [0.78519918, 2.25476522, 4.16710899, 11.83677242, 2.46256654, 5.49154145, 4.77488592, 10.82382343, 2.55091774, 0.0, 11.84912271, 12.07478086, 10.29278088, 6.0790813, 13.94477075, 14.95145432, 4.69860935, 4.35512136, 7.68864491, 7.1942034, 2.59516155, 3.84654483, 6.20912235, 6.41300434, 15.6272628, 5.62226734], [12.41331877, 9.73701553, 16.01513254, 1.21847207, 12.99438637, 6.46790425, 13.11344175, 2.09441163, 10.8229774, 11.84912271, 0.0, 1.90003438, 1.56357422, 17.13358255, 4.02199744, 5.66167899, 7.29470888, 15.5678518, 4.89030374, 4.69686133, 14.4410741, 8.60373002, 12.71754856, 14.28325341, 5.25445843, 16.03227189], [12.71590282, 10.12917852, 16.21524974, 3.11832915, 12.87478936, 6.58377272, 12.62301977, 3.84432566, 10.67314007, 12.07478086, 1.90003438, 0.0, 2.6118074, 17.70147714, 2.3518661, 3.93205427, 7.83263402, 16.06600072, 5.99873107, 5.31296985, 14.65862467, 8.50059261, 11.98375994, 13.61875396, 3.86709812, 16.71557835], [10.85152542, 8.17368452, 14.45748061, 1.83175951, 11.49295696, 4.9515114, 11.7338606, 1.64695079, 9.33868012, 10.29278088, 1.56357422, 2.6118074, 0.0, 15.58271487, 4.96315262, 6.53297558, 5.73233014, 14.00948474, 3.45819553, 3.13352906, 12.88307769, 7.12195452, 11.45710668, 12.9758221, 6.42125175, 14.49783156], [5.30531019, 7.60734617, 3.5593164, 16.84611761, 7.07799438, 11.23322217, 9.50479511, 15.6963289, 8.46980964, 6.0790813, 17.13358255, 17.70147714, 15.58271487, 0.0, 19.74369024, 20.86788685, 9.8867652, 1.72409828, 12.40949609, 12.47385352, 4.15383123, 9.89548123, 11.26718792, 10.69383476, 21.41887369, 1.50735676], [14.62988656, 12.13698276, 18.01879751, 5.18439858, 14.46175496, 8.51841325, 13.85853526, 6.09367811, 12.27957352, 13.94477075, 4.02199744, 2.3518661, 4.96315262, 19.74369024, 0.0, 1.64194246, 9.96019735, 18.07538198, 8.33218791, 7.53462099, 16.49043977, 10.20833345, 12.95956846, 14.63838209, 1.68260264, 18.82844936], [15.66729075, 13.27478887, 18.94081739, 6.82614445, 15.24009862, 9.64893718, 14.38249205, 7.72561094, 13.1003455, 14.95145432, 5.66167899, 3.93205427, 6.53297558, 20.86788685, 1.64194246, 0.0, 11.22266398, 19.17578297, 9.82799986, 8.90862517, 17.44571855, 11.13886346, 13.31099511, 14.99462349, 1.51257476, 20.02425449], [5.1738933, 2.47957666, 8.81211911, 7.17018836, 6.40054883, 1.89886861, 7.54307763, 6.12954017, 4.58428986, 4.69860935, 7.29470888, 7.83263402, 5.73233014, 9.8867652, 9.96019735, 11.22266398, 0.0, 8.28568492, 3.02503932, 2.60056235, 7.24217467, 2.90742183, 8.06758082, 9.13380464, 11.61420387, 8.88477454], [3.58142053, 5.94352464, 2.56523183, 15.33762043, 5.49462281, 9.55770251, 7.98929774, 14.21040127, 6.76531945, 4.35512136, 15.5678518, 16.06600072, 14.00948474, 1.72409828, 18.07538198, 19.17578297, 8.28568492, 0.0, 10.92506501, 10.88435265, 2.61067234, 8.17412445, 9.72321638, 9.28935165, 19.75367659, 1.74174469], [8.10023862, 5.43828894, 11.73700783, 4.44343542, 9.40801035, 3.47662582, 10.34049008, 3.28797438, 7.4908784, 7.68864491, 4.89030374, 5.99873107, 3.45819553, 12.40949609, 8.33218791, 9.82799986, 3.02503932, 10.92506501, 0.0, 1.58089135, 10.18172332, 5.50726129, 10.59249795, 11.84585271, 9.85510696, 11.23092156], [7.72959113, 5.04438059, 11.34952072, 4.65384074, 8.57909014, 2.17386813, 9.1927582, 3.72233099, 6.51112833, 7.1942034, 4.69686133, 5.31296985, 3.13352906, 12.47385352, 7.53462099, 8.90862517, 2.60056235, 10.88435265, 1.58089135, 0.0, 9.77500048, 4.369047, 9.26792185, 10.61550814, 9.1523327, 11.43006949], [2.08169749, 4.76417279, 1.57468355, 14.40544845, 2.92421481, 8.07913676, 5.38282355, 13.37009259, 4.51501243, 2.59516155, 14.4410741, 14.65862467, 12.88307769, 4.15383123, 16.49043977, 17.44571855, 7.24217467, 2.61067234, 10.18172332, 9.77500048, 0.0, 6.3072182, 7.1281481, 6.68967693, 18.17146896, 4.30327841], [4.59862711, 2.78090837, 7.8158263, 8.8601632, 4.40130551, 2.21692624, 4.84396878, 8.05330675, 2.21931093, 3.84654483, 8.60373002, 8.50059261, 7.12195452, 9.89548123, 10.20833345, 11.13886346, 2.90742183, 8.17412445, 5.50726129, 4.369047, 6.3072182, 0.0, 5.17553967, 6.34196185, 11.88525069, 9.28096237], [6.84565681, 7.02480049, 7.95886858, 13.28726193, 4.23596374, 7.13140097, 1.78281258, 12.71380328, 3.98174339, 6.20912235, 12.71754856, 11.98375994, 11.45710668, 11.26718792, 12.95956846, 13.31099511, 8.06758082, 9.72321638, 10.59249795, 9.26792185, 7.1281481, 5.17553967, 0.0, 1.68452525, 14.47827637, 11.3473849], [6.91266051, 7.70310563, 7.21101714, 14.79876271, 4.06829594, 8.44665983, 1.68452975, 14.16409044, 4.67500764, 6.41300434, 14.28325341, 13.61875396, 12.9758221, 10.69383476, 14.63838209, 14.99462349, 9.13380464, 9.28935165, 11.84585271, 10.61550814, 6.68967693, 6.34196185, 1.68452525, 0.0, 16.16189478, 10.99268415], [16.3123257, 13.81340365, 19.69790598, 6.2995567, 16.11947556, 10.19802227, 15.44610294, 7.34540194, 13.94327444, 15.6272628, 5.25445843, 3.86709812, 6.42125175, 21.41887369, 1.68260264, 1.51257476, 11.61420387, 19.75367659, 9.85510696, 9.1523327, 18.17146896, 11.88525069, 14.47827637, 16.16189478, 0.0, 20.49131732], [4.8405079, 6.77107657, 4.25692634, 15.67399502, 7.11223997, 10.37950982, 9.64217133, 14.50372216, 8.14568881, 5.62226734, 16.03227189, 16.71557835, 14.49783156, 1.50735676, 18.82844936, 20.02425449, 8.88477454, 1.74174469, 11.23092156, 11.43006949, 4.30327841, 9.28096237, 11.3473849, 10.99268415, 20.49131732, 0.0]]
            
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
            for j in range(T):
                temp.append(0)
            b.append(temp)
        for i in range(n):
            for j in range(T):
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
        p_prime = []
        for i in range(n):
            temp = []
            for j in range(T):
                temp.append(0)
            p_prime.append(temp)
        for i in range(n):
            for j in range(T):
                p_prime[i][j] = m.addVar(vtype=GRB.BINARY, name="p_prime,%s" % str(i+1) + "," + str(j+1))                
                
        # Con esta garantizo que p_prime copie lo de p (los 1s)
        for j in range(T):
            for i in range(n):
                m.addConstr(p_prime[i][j] >= p[i][j])
                
        # esta garantiza que haya sume a 1 las columnas de p_prime
        for j in range(T):
            sum_ = 0
            for i in range(n):
                sum_ = sum_ + p_prime[i][j]
            m.addConstr(sum_ == 1)
            
        # Con esta, las columnas que suman cero en p, se convierten en p_prime = p_prime anterior
        for j in range(T):
            for i in range(n):
                sum_ = 0
                for k in range(n):
                    sum_ = sum_ + p[k][j]
                if j == 0:
                    m.addConstr(p_prime[i][j] >= 0) # Deshabilitada, porque siempre defiende a alguien al inicio
                else:
                    m.addConstr(p_prime[i][j] >= p_prime[i][j-1] * (1 - sum_))
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
            for j in range(T):
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
       
        
        for j in range(T): # para evitar que se defiendan más de 1
            sum_ = 0
            for i in range(n):
                sum_ = sum_ + p[i][j]
            m.addConstr(sum_ <= 1)
     
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
                        sum_1 = sum_1 + p_prime[k][j] * mov[k][i]
                    dist_nodos = dist_nodos + sum_1 * d0[i]
                else:
                    for k in range(n):
                        sum_1 = sum_1 + p_prime[k][j] * mov[k][i]
                    dist_nodos = dist_nodos + sum_1 * p_prime[i][j-1]
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
        sum_ = 0
        for i in range(n):
            sum_ = sum_ + p[i][0]
        m.addConstr(sum_ == 1)
        # Y se deben defender siempre en cada ronda (consecutivamente)
        for j in range(1,T):
            sum_1 = 0
            sum_2 = 0
            for i in range(n):
                sum_1 = sum_1 + p[i][j-1]
                sum_2 = sum_2 + p[i][j]
            m.addConstr(sum_2 <= sum_1)
 
        # Manejar T como entrada del problema y checar la solución
        # para ver si es consistente y minimiza los quemados
 
        #for j in range(T):
        for j in range(Td): # Hay que encontrar este valor en automático
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
                
        #m.addConstr(d[8][1] == 1)
                                        
        #---------------------------- OBJECTIVE FUNCTION --------------------------------------------------        
        b_transpose = np.array(b).T.tolist()
        #m.setObjective(sum(b_transpose[T-1]), GRB.MINIMIZE)#-----------------------------------------------(1)
        d_transpose = np.array(d).T.tolist()
        m.setObjective(sum(d_transpose[T-1]), GRB.MAXIMIZE)#-----------------------------------------------(1)
        #---------------------------- OPTIMIZATION -------------------------------------------------------
        
        #m.computeIIS()
        #m.write("model.lp")
        
        
        m.optimize()
        runtime = m.Runtime
        print("The run time is %f" % runtime)
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
            #b[j][k] = b_out[i]
            k =  k + 1
            if (i+1) % T == 0:
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
        feasible = True
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
            best_defense      = defense_out[:Td]
            best_defense_time = time_out[:Td]
            
                
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
    global instance, feasible, best_defense, best_defense_time
    folder_dataset = 'C:/Users/jgd/Documents/FireFighter/'
    dataset = [
        #['lattice4x4.mtx', 16, 24, [15]]
        ['tree.mtx', 26, 24, [18]]
        #['lattice2x2.mtx', 4, 4, [3]]
        #['path4.mtx', 4, 3, [3]]
        #['test_lattice2x2.mtx', 5, 5, [2]]
        #['test2_lattice2x2.mtx', 6, 6, [2]]
        #['test3_lattice2x2.mtx', 8, 8, [2]]
        #['path6.mtx', 6, 5, [5]]
        ]
    T = 20
    T2 = 5
    best_defense = []
    best_defense_time = []
    for i in range(len(dataset)):
        print("--------------------------------------------------------------")
        instance = dataset[i][0]
        print("instance: " + instance)
        main(dataset[i][1], dataset[i][2], folder_dataset + dataset[i][0], dataset[i][3])
        upper = 10
        lower = 1
        while lower <= upper:
            mid = math.floor((upper + lower) /2)
            print("MID: " + str(mid))
            run(T, T2, mid)
            #total_runtime = total_runtime + runtime
            if feasible:
                upper = mid - 1
            else:
                lower = mid + 1
        print(best_defense)
        print(best_defense_time)
