--#########################################################################
--##                                                                     ##
--##  DESLOCA: Uma ferramenta gráfico-interativo-didática para o ensino  ##
--##  do Método dos Deslocamentos em estruturas bidimensionais           ##
--##                                                                     ##
--##  Autor: Bruno Valadão Cunha - 2018                                  ##
--##                                                                     ##
--##  Desenvolvido em Lua https://www.lua.org/ e utilizando as           ##
--##  bibliotecas wxLua http://wxlua.sourceforge.net/ e                  ##
--##  e LuaMatrix https://github.com/davidm/lua-matrix                   ##
--##                                                                     ##
--##  Trabalho de Conclusão de Curso em Engenharia Civil - UFRGS         ##
--##                                                                     ##
--##  Repositório: https://github.com/brunovalads/DESLOCA                ##
--##                                                                     ##
--#########################################################################


--############################################################################################################################################################
-- DESLOCA
--############################################################################################################################################################

-- Import external module
local matrix = require "matrix"

-- DESLOCA's global functions and constants
local DESLOCA = {}

-- DESLOCA's version number
DESLOCA.version = "0.28"

-- Structure input, each entry is a bar, following the pattern xi, yi, xf, yf, E, A, I, dist {q, dist, is_global, global_dir} (in cm, kN/cm^2, cm^2, cm^4, kN/cm)
DESLOCA.bar = { -- TODO: see if this will be the deafult structure or the user will start from scratch
  --    xi,  yi,  xf,  yf,     E,  A,     I
  [1] = {xi = 0, yi = 0, xf = 0, yf = 400, E = 20000, A = 60, I = 8000, global_node = {}, dist = {q = -0.2, is_global = false, global_dir = "x"}, -- bar 1
         i = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = false, free_y = false, free_z = false},
         f = {Fx = 50, Fy = 0, Mz = 3000, is_global = true, free_x = true, free_y = true, free_z = true}},
  
  [2] = {xi = 0, yi = 400, xf = 600, yf = 400, E = 20000, A = 80, I = 16000, global_node = {}, dist = {q = -0.1, is_global = false}, -- bar 2
         i = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = true, free_y = true, free_z = true},
         f = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = false, free_y = false, free_z = false}},
}

-- Make a copy of the default structure
DESLOCA.default = {
  [1] = {xi = 0, yi = 0, xf = 0, yf = 400, E = 20000, A = 60, I = 8000, global_node = {}, dist = {q = -0.2, is_global = false, global_dir = "x"}, -- bar 1
         i = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = false, free_y = false, free_z = false},
         f = {Fx = 50, Fy = 0, Mz = 3000, is_global = true, free_x = true, free_y = true, free_z = true}},
  [2] = {xi = 0, yi = 400, xf = 600, yf = 400, E = 20000, A = 80, I = 16000, global_node = {}, dist = {q = -0.1, is_global = false}, -- bar 2
         i = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = true, free_y = true, free_z = true},
         f = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = false, free_y = false, free_z = false}},
}

-- Function that returns the size and angle of the bar
function DESLOCA.dimension(bar)

  local delta_x = bar.xf - bar.xi
  local delta_y = bar.yf - bar.yi

  --if delta_x == 0 and delta_y == 0 then error("Alguma barra é um ponto!") end -- TODO: replace error for message
  
  local L = math.sqrt(delta_x^2 + delta_y^2)
  
  local cosine = delta_x/L
  local sine = delta_y/L
  
  return L, cosine, sine
 
end

-- Function that calculates and returns the stiffness matrix of the bar in local and global coordinates and its rotation matrixes
function DESLOCA.get_bar_matrix(bar)
  
  -- Import bar properties
  local E, A, I = bar.E, bar.A, bar.I
  local L, C, S = DESLOCA.dimension(bar)
  
  -- Stiffness coefficients
  local k1 = E*A/L
  local k2 = 12*E*I/(L^3)
  local k3 = 6*E*I/(L^2)
  local k4 = 4*E*I/L
  local k5 = 2*E*I/L
  
  -- Construction of the stiffness matrix in local coordinates
  local stiff_matrix_local = {
    { k1,  0,  0,-k1,  0,  0},
    {  0, k2, k3,  0,-k2, k3},
    {  0, k3, k4,  0,-k3, k5},
    {-k1,  0,  0, k1,  0,  0},
    {  0,-k2,-k3,  0, k2,-k3},
    {  0, k3, k5,  0,-k3, k4}
  }
  
  -- Construction of the rotation matrix -- TODO: not used yet, but interesting for didactic purposes
  local rotation_matrix = {
    { C, S, 0, 0, 0, 0},
    {-S, C, 0, 0, 0, 0},
    { 0, 0, 1, 0, 0, 0},
    { 0, 0, 0, C, S, 0},
    { 0, 0, 0,-S, C, 0},
    { 0, 0, 0, 0, 0, 1}
  }
  
  -- Construction of the transposed rotation matrix -- TODO: not used yet, but interesting for didactic purposes
  local rotation_matrix_t = {
    { C,-S, 0, 0, 0, 0},
    { S, C, 0, 0, 0, 0},
    { 0, 0, 1, 0, 0, 0},
    { 0, 0, 0, C,-S, 0},
    { 0, 0, 0, S, C, 0},
    { 0, 0, 0, 0, 0, 1}
  }
  
  -- Construction of the stiffness matrix in global coordinates
  local stiff_matrix_global = {
    {  k1*C^2+k2*S^2,    (k1-k2)*C*S, -k3*S, -k1*C^2-k2*S^2,    (k2-k1)*C*S, -k3*S},
    {    (k1-k2)*C*S,  k1*S^2+k2*C^2,  k3*C,    (k2-k1)*C*S, -k1*S^2-k2*C^2,  k3*C},
    {          -k3*S,           k3*C,    k4,           k3*S,          -k3*C,    k5},
    { -k1*C^2-k2*S^2,    (k2-k1)*C*S,  k3*S,  k1*C^2+k2*S^2,    (k1-k2)*C*S,  k3*S},
    {    (k2-k1)*C*S, -k1*S^2-k2*C^2, -k3*C,    (k1-k2)*C*S,  k1*S^2+k2*C^2, -k3*C},
    {          -k3*S,           k3*C,    k5,           k3*S,          -k3*C,    k4}
  }
  
  return stiff_matrix_local, stiff_matrix_global, rotation_matrix, rotation_matrix_t
end

-- Function that calculates and return the load vector of the bar in local and global coordinates
function DESLOCA.get_bar_load_vector(bar)
  
  -- Import bar properties
  local L, C, S = DESLOCA.dimension(bar) 
  local q_local, q_global
  local is_global = bar.dist.is_global
  if is_global then
    q_global = bar.dist.q
    if bar.dist.global_dir == "x" then
      q_local = q_global*S
    elseif bar.dist.global_dir == "y" then
      q_local = q_global*C
    end
  else
    q_local = bar.dist.q
    q_global = q_local*C -- this one is never used here    
  end
  
  -- Load coefficients
  local Fxi = 0
  local Fyi = q_local*L/2
  local Mzi = q_local*(L^2)/12
  local Fxf = 0
  local Fyf = q_local*L/2
  local Mzf = -q_local*(L^2)/12
  
  -- Construction of the load vector in local coordinates
  local f_local = {Fxi, Fyi, Mzi, Fxf, Fyf, Mzf}
  
  -- Construction of the transposed rotation matrix
  local Rt = {
    { C,-S, 0, 0, 0, 0},
    { S, C, 0, 0, 0, 0},
    { 0, 0, 1, 0, 0, 0},
    { 0, 0, 0, C,-S, 0},
    { 0, 0, 0, S, C, 0},
    { 0, 0, 0, 0, 0, 1}
  }
  
  -- Construction of the load vector in global coordinates
  local f_global = {
    Rt[1][1]*f_local[1] + Rt[1][2]*f_local[2] + Rt[1][3]*f_local[3] + Rt[1][4]*f_local[4] + Rt[1][5]*f_local[5] + Rt[1][6]*f_local[6],
    Rt[2][1]*f_local[1] + Rt[2][2]*f_local[2] + Rt[2][3]*f_local[3] + Rt[2][4]*f_local[4] + Rt[2][5]*f_local[5] + Rt[2][6]*f_local[6],
    Rt[3][1]*f_local[1] + Rt[3][2]*f_local[2] + Rt[3][3]*f_local[3] + Rt[3][4]*f_local[4] + Rt[3][5]*f_local[5] + Rt[3][6]*f_local[6],
    Rt[4][1]*f_local[1] + Rt[4][2]*f_local[2] + Rt[4][3]*f_local[3] + Rt[4][4]*f_local[4] + Rt[4][5]*f_local[5] + Rt[4][6]*f_local[6],
    Rt[5][1]*f_local[1] + Rt[5][2]*f_local[2] + Rt[5][3]*f_local[3] + Rt[5][4]*f_local[4] + Rt[5][5]*f_local[5] + Rt[5][6]*f_local[6],
    Rt[6][1]*f_local[1] + Rt[6][2]*f_local[2] + Rt[6][3]*f_local[3] + Rt[6][4]*f_local[4] + Rt[6][5]*f_local[5] + Rt[6][6]*f_local[6]
  }
  
  --[[ DEBUG REMOVE/TEST
  print("----------------------------------------------")
  print(string.format("\nis_global: %s\nq_local: %.3f, q_global: %.3f", is_global and "true" or "false", q_local, q_global))
  
  print("\nRt = ")
  for i = 1, 6 do
    print(Rt[i][1], Rt[i][2], Rt[i][3], Rt[i][4], Rt[i][5], Rt[i][6])
  end
  
  print("\nf_local = ")
  print(f_local[1], f_local[2], f_local[3], f_local[4], f_local[5], f_local[6])
  
  print("\nf_global = ")
  print(f_global[1], f_global[2], f_global[3], f_global[4], f_global[5], f_global[6])]]
  
  return f_local, f_global  
end

-- Function that calculates and returns the stiffness matrix of the structure
function DESLOCA.get_structure_stiff_matrix()

  -- Ensure bar global matrixes are calculated
  if DESLOCA.bar[1].stiff_matrix_global == nil then
    DESLOCA.all_matrix_update() -- TODO: all_matrix_update() has get_structure_stiff_matrix(), so it calculates twice if this call happens
  end
  
  -- Import bar info and divide it in supercells
  local bar_global_matrixes = {}
  for i = 1, #DESLOCA.bar do
    bar_global_matrixes[i] = {}
    bar_global_matrixes[i].Kii, bar_global_matrixes[i].Kif, bar_global_matrixes[i].Kfi, bar_global_matrixes[i].Kff = {}, {}, {}, {}
    
    -- Kii supercell
    for j = 1, 3 do
      bar_global_matrixes[i].Kii[j] = {}
      for k = 1, 3 do
        bar_global_matrixes[i].Kii[j][k] = DESLOCA.bar[i].stiff_matrix_global[j][k]
      end
    end
    -- Kif supercell
    for j = 1, 3 do
      bar_global_matrixes[i].Kif[j] = {}
      for k = 4, 6 do
        bar_global_matrixes[i].Kif[j][k-3] = DESLOCA.bar[i].stiff_matrix_global[j][k]
      end
    end
    -- Kfi supercell
    for j = 4, 6 do
      bar_global_matrixes[i].Kfi[j-3] = {}
      for k = 1, 3 do
        bar_global_matrixes[i].Kfi[j-3][k] = DESLOCA.bar[i].stiff_matrix_global[j][k]
      end
    end
    -- Kff supercell
    for j = 4, 6 do
      bar_global_matrixes[i].Kff[j-3] = {}
      for k = 4, 6 do
        bar_global_matrixes[i].Kff[j-3][k-3] = DESLOCA.bar[i].stiff_matrix_global[j][k]
      end
    end
    
  end
  
  -- Init the structure stiffness matrixes
  local K = {}
  local K_simp = {}
  for i = 1, 3*#DESLOCA.node do
    K[i] = {}
    K_simp[i] = {}
    for j = 1, 3*#DESLOCA.node do
      K[i][j] = 0
      K_simp[i][j] = 0
    end
  end
  
  -- Fill both matrixes (K_simp will be simplified later)
  for i = 1, #DESLOCA.bar do
    local node_i, node_f = DESLOCA.bar[i].global_node.i, DESLOCA.bar[i].global_node.f
    
    -- ii
    for j = 3*node_i - 2, 3*node_i do
      for k = 3*node_i - 2, 3*node_i do
        K[j][k] = K[j][k] + bar_global_matrixes[i].Kii[j-3*(node_i-1)][k-3*(node_i-1)]
        K_simp[j][k] = K[j][k] -- this is done to copy K into K_simp, since there's no table copy in Lua and matrix.copy (module) is slow
      end
    end
    -- if
    for j = 3*node_i - 2, 3*node_i do
      for k = 3*node_f - 2, 3*node_f do
        K[j][k] = K[j][k] + bar_global_matrixes[i].Kif[j-3*(node_i-1)][k-3*(node_f-1)]
        K_simp[j][k] = K[j][k]
      end
    end
    -- fi
    for j = 3*node_f - 2, 3*node_f do
      for k = 3*node_i - 2, 3*node_i do
        K[j][k] = K[j][k] + bar_global_matrixes[i].Kfi[j-3*(node_f-1)][k-3*(node_i-1)]
        K_simp[j][k] = K[j][k]
      end
    end
    -- ff
    for j = 3*node_f - 2, 3*node_f do
      for k = 3*node_f - 2, 3*node_f do
        K[j][k] = K[j][k] + bar_global_matrixes[i].Kff[j-3*(node_f-1)][k-3*(node_f-1)]
        K_simp[j][k] = K[j][k]
      end
    end
    
  end
  
  -- Simplify the structure simplified stiffness matrix
  for j = 1, 3*#DESLOCA.node do
    for k = 1, 3*#DESLOCA.node do
    
      if (not DESLOCA.node[math.ceil(k/3)].free_x and (k-1)%3 + 1 == 1) or (not DESLOCA.node[math.ceil(j/3)].free_x and (j-1)%3 + 1 == 1) then
        
        K_simp[j][k] = j == k and 1 or 0 -- 1 if diagonal, 0 if not
        
      end
      if (not DESLOCA.node[math.ceil(k/3)].free_y and (k-1)%3 + 1 == 2) or (not DESLOCA.node[math.ceil(j/3)].free_y and (j-1)%3 + 1 == 2) then
        
        K_simp[j][k] = j == k and 1 or 0 -- 1 if diagonal, 0 if not
        
      end
      if (not DESLOCA.node[math.ceil(k/3)].free_z and (k-1)%3 + 1 == 3) or (not DESLOCA.node[math.ceil(j/3)].free_z and (j-1)%3 + 1 == 3) then
        
        K_simp[j][k] = j == k and 1 or 0 -- 1 if diagonal, 0 if not
        
      end
      
    end
  end
  
  return K, K_simp
end

-- Function that calculates and return the load vector of the structure in global coordinates
function DESLOCA.get_structure_load_vector()
  
  -- Ensure bar global load vectors are calculated
  if DESLOCA.bar[1].load_vector_global == nil then
    DESLOCA.load_vector_update() -- TODO: load_vector_update() has get_structure_load_vector(), so it calculates twice if this call happens
  end
  
  -- Import bar info and divide it in sub-vectors
  local bar_global_vectors = {}
  for i = 1, #DESLOCA.bar do
    bar_global_vectors[i] = {}
    bar_global_vectors[i].Fi, bar_global_vectors[i].Ff = {}, {}
    
    -- Fi sub-vector
    for j = 1, 3 do
      bar_global_vectors[i].Fi[j] = DESLOCA.bar[i].load_vector_global[j]
    end
    -- Ff sub-vector
    for j = 4, 6 do
      bar_global_vectors[i].Ff[j-3] = DESLOCA.bar[i].load_vector_global[j]
    end
    
  end
  
  -- Init the structure load vectors
  local F, F_simp = {}, {} -- (F is the result vector, F_simp is the vector with boundary conditions)
  local F_dist, F_nodal = {}, {} -- (F_dist and F_nodal are the portions due distributed load and nodal load respectively)
  for i = 1, 3*#DESLOCA.node do
    F[i] = 0
    F_simp[i] = 0
    F_dist[i] = 0
    F_nodal[i] = 0
  end
  
  -- Fill load vectors (F_simp will be simplified later)
  for i = 1, #DESLOCA.bar do
    local node_i, node_f = DESLOCA.bar[i].global_node.i, DESLOCA.bar[i].global_node.f
    
    -- i
    for j = 3*node_i - 2, 3*node_i do
      F[j] = F[j] + bar_global_vectors[i].Fi[j-3*(node_i-1)]
      F_simp[j] = F[j]
      F_dist[j] = F[j]
    end
    -- f
    for j = 3*node_f - 2, 3*node_f do
      F[j] = F[j] + bar_global_vectors[i].Ff[j-3*(node_f-1)]
      F_simp[j] = F[j]
      F_dist[j] = F[j]
    end
    
  end
  
  -- Add the nodal loads in both load vectors (F_simp will be simplified later)
  for j = 1, 3*#DESLOCA.node do
    
    if (j-1)%3 + 1 == 1 then
      F[j] = F[j] + DESLOCA.node[math.ceil(j/3)].Fx
      F_simp[j] = F[j]
      F_nodal[j] = F_nodal[j] + DESLOCA.node[math.ceil(j/3)].Fx
    elseif (j-1)%3 + 1 == 2 then
      F[j] = F[j] + DESLOCA.node[math.ceil(j/3)].Fy
      F_simp[j] = F[j]
      F_nodal[j] = F_nodal[j] + DESLOCA.node[math.ceil(j/3)].Fy
    elseif (j-1)%3 + 1 == 3 then
      F[j] = F[j] + DESLOCA.node[math.ceil(j/3)].Mz
      F_simp[j] = F[j]
      F_nodal[j] = F_nodal[j] + DESLOCA.node[math.ceil(j/3)].Mz
    end
    
  end
  
  -- Simplify the structure simplified load vector
  for j = 1, 3*#DESLOCA.node do
    
    if not DESLOCA.node[math.ceil(j/3)].free_x and (j-1)%3 + 1 == 1 then
      F_simp[j] = 0
    end
    if not DESLOCA.node[math.ceil(j/3)].free_y and (j-1)%3 + 1 == 2 then
      F_simp[j] = 0
    end
    if not DESLOCA.node[math.ceil(j/3)].free_z and (j-1)%3 + 1 == 3 then
      F_simp[j] = 0
    end
      
  end
  
  --[[ DEBUG
  print(string.format("\nStructure load vector (%d memebers)", #F))
  for i = 1, #F do
    print(F[i])
  end]]
  
  return F, F_simp, F_dist, F_nodal
end

-- Update all the matrixes
function DESLOCA.all_matrix_update() -- TODO: rename if Structure matrix remains inside here
  -- Bar matrixes
  for i = 1, #DESLOCA.bar do
    DESLOCA.bar[i].stiff_matrix_local, DESLOCA.bar[i].stiff_matrix_global, DESLOCA.bar[i].rotation_matrix, DESLOCA.bar[i].rotation_matrix_t = DESLOCA.get_bar_matrix(DESLOCA.bar[i])
  end
  
  -- Structure matrixes
  DESLOCA.structure_stiff_matrix, DESLOCA.structure_stiff_matrix_simp = DESLOCA.get_structure_stiff_matrix()
end

-- Update all the load vectors
function DESLOCA.load_vector_update()
  -- Bar vectors
  for i = 1, #DESLOCA.bar do
    DESLOCA.bar[i].load_vector_local, DESLOCA.bar[i].load_vector_global = DESLOCA.get_bar_load_vector(DESLOCA.bar[i])
  end
  
  -- Structure vectors
  DESLOCA.structure_load_vector, DESLOCA.structure_load_vector_simp, DESLOCA.structure_load_vector_dist_portion, DESLOCA.structure_load_vector_nodal_portion = DESLOCA.get_structure_load_vector()
end

-- Function that scans the structure enumerating the nodes
function DESLOCA.eval_nodes()
  DESLOCA.node = {} -- need to reset every time

  local structure_size = #DESLOCA.bar
  if structure_size < 1 then error("Insira pelo menos 1 barra na estrutura!") end -- TODO: replace error for message (probably this will never happen anyways)
  
  -- Initial nodes, from the 1st bar
  if math.round(DESLOCA.bar[1].xi) == math.round(DESLOCA.bar[1].xf) and math.round(DESLOCA.bar[1].yi) == math.round(DESLOCA.bar[1].yf) then error("A barra 1 é um ponto!") end -- TODO: replace error for message (probably this will never happen anyways)
  DESLOCA.node[1] = {x = DESLOCA.bar[1].xi, y = DESLOCA.bar[1].yi, bars = {1}}
  DESLOCA.node[2] = {x = DESLOCA.bar[1].xf, y = DESLOCA.bar[1].yf, bars = {1}}
  
  -- Node enumerating on the bars 2+
  if structure_size > 1 then
  
    -- Nodes before the equals overlay
    local pre_nodes = {}
    pre_nodes[1] = DESLOCA.node[1]
    pre_nodes[2] = DESLOCA.node[2]
    
    for i = 2, structure_size do
      
      if math.round(DESLOCA.bar[i].xi) == math.round(DESLOCA.bar[i].xf) and math.round(DESLOCA.bar[i].yi) == math.round(DESLOCA.bar[i].yf) then error("A barra " .. i .. " é um ponto!") end -- TODO: replace error for message (probably this will never happen anyways)
      
      table.insert(pre_nodes, {x = DESLOCA.bar[i].xi, y = DESLOCA.bar[i].yi, bars = {i}})
      table.insert(pre_nodes, {x = DESLOCA.bar[i].xf, y = DESLOCA.bar[i].yf, bars = {i}})
      
    end
    
    -- Node overlay
    for i = 3, #pre_nodes do
      
      local overlay = false
      
      for j = 1, i - 1 do 
      
        if math.round(pre_nodes[i].x) == math.round(pre_nodes[j].x) and math.round(pre_nodes[i].y) == math.round(pre_nodes[j].y) then -- there is node overlay, excluding decimals of centimeters
          
          overlay = true
          table.insert(pre_nodes[i].bars, math.floor((j+1)/2))
          table.insert(pre_nodes[j].bars, math.floor((i+1)/2))
          break
          
        else -- is a new node
          
          overlay = false
          
        end
        
      end
      
      if not overlay then
      
        table.insert(DESLOCA.node, {x = pre_nodes[i].x, y = pre_nodes[i].y, bars = pre_nodes[i].bars})
      
      end
      
    end
    
  end
  
  -- Indentify if nodes are the initial or final  
  for j = 1, #DESLOCA.node do
    DESLOCA.node[j].bar_nodes = {}
  end
  for i = 1, structure_size do
    
    for j = 1, #DESLOCA.node do
      
      if math.round(DESLOCA.bar[i].xi) == math.round(DESLOCA.node[j].x) and math.round(DESLOCA.bar[i].yi) == math.round(DESLOCA.node[j].y) then
        table.insert(DESLOCA.node[j].bar_nodes, "i")
        DESLOCA.bar[i].global_node.i = j
      elseif math.round(DESLOCA.bar[i].xf) == math.round(DESLOCA.node[j].x) and math.round(DESLOCA.bar[i].yf) == math.round(DESLOCA.node[j].y) then
        table.insert(DESLOCA.node[j].bar_nodes, "f")
        DESLOCA.bar[i].global_node.f = j
      end
    end
  end
  
  -- Import restraints
  for j = 1, #DESLOCA.node do
    DESLOCA.node[j].free_x = true
    DESLOCA.node[j].free_y = true
    DESLOCA.node[j].free_z = true
    for k = 1, #DESLOCA.node[j].bars do
      if DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].free_x == false then DESLOCA.node[j].free_x = false end
      if DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].free_y == false then DESLOCA.node[j].free_y = false end
      if DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].free_z == false then DESLOCA.node[j].free_z = false end
    end
  end
  
  -- Import nodal loads
  for j = 1, #DESLOCA.node do
    DESLOCA.node[j].Fx = 0
    DESLOCA.node[j].Fy = 0
    DESLOCA.node[j].Mz = 0
    for k = 1, #DESLOCA.node[j].bars do
      
      if DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].is_global then -- just add the loads
        DESLOCA.node[j].Fx = DESLOCA.node[j].Fx + DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].Fx
        DESLOCA.node[j].Fy = DESLOCA.node[j].Fy + DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].Fy
        DESLOCA.node[j].Mz = DESLOCA.node[j].Mz + DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].Mz
      else -- manipulate the loads to properly add them
        -- Import bar properties
        local L, C, S = DESLOCA.dimension(DESLOCA.bar[DESLOCA.node[j].bars[k]])
        local angle_rad = Get_angle(S, C)
        local bar_Fx = DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].Fx -- relative to the node ("i" or "f") stored in .bar_nodes[k]
        local bar_Fy = DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].Fy -- relative to the node ("i" or "f") stored in .bar_nodes[k]
        local bar_Mz = DESLOCA.bar[DESLOCA.node[j].bars[k]][DESLOCA.node[j].bar_nodes[k]].Mz -- relative to the node ("i" or "f") stored in .bar_nodes[k]
        
        DESLOCA.node[j].Fx = DESLOCA.node[j].Fx + bar_Fx*C + bar_Fy*math.cos(angle_rad + math.pi/2)
        DESLOCA.node[j].Fy = DESLOCA.node[j].Fy + bar_Fx*S + bar_Fy*math.sin(angle_rad + math.pi/2)
        DESLOCA.node[j].Mz = DESLOCA.node[j].Mz + bar_Mz -- this one doesn't matter if global or local
      end
    end
  end
  
  
  --[[
  -- DEBUG -- REMOVE/TEST
  print("\nListar barras:")
  for i = 1, #DESLOCA.bar do
    
    local bar_str = string.format("Barra %d = {xi = %d, yi = %d, xf = %d, yf = %d, global_node = {%d, %d}}", i, DESLOCA.bar[i].xi, DESLOCA.bar[i].yi, DESLOCA.bar[i].xf, DESLOCA.bar[i].xf, DESLOCA.bar[i].global_node.i, DESLOCA.bar[i].global_node.f)
    print(bar_str)
  end
  --
  print("\nListar nós:")
  for j = 1, #DESLOCA.node do
  
    local nos_barras_str, bar_nodes_str = "", ""
    for i = 1, #DESLOCA.node[j].bars do
      nos_barras_str = nos_barras_str .. DESLOCA.node[j].bars[i]
      bar_nodes_str = bar_nodes_str .. DESLOCA.node[j].bar_nodes[i]
      if i ~= #DESLOCA.node[j].bars then nos_barras_str = nos_barras_str .. ", " ; bar_nodes_str = bar_nodes_str .. ", " end
    end
    
    local nos_str = string.format("Nó %d = {x = %d, y = %d, barras = {%s}, bar_nodes = {%s}, free_x = %s, free_y = %s, free_z = %s, Fx = %.1f, Fy = %.1f, Mz = %.1f}, ",
    j, DESLOCA.node[j].x, DESLOCA.node[j].y, nos_barras_str, bar_nodes_str, DESLOCA.node[j].free_x and "true" or "false", DESLOCA.node[j].free_y and "true" or "false", DESLOCA.node[j].free_z and "true" or "false",
    DESLOCA.node[j].Fx, DESLOCA.node[j].Fy, DESLOCA.node[j].Mz)
    print(nos_str)
  end]]
  
end

-- Function that prints a matrix in the console, just sed for debugging
local function print_mat(mat)
  print(string.format("\n%dx%d matrix:", #mat, #mat))
  local wolfram_str = "{" -- to check results in Wolfram Alpha
  for i = 1, #mat do
    if i == 1 then wolfram_str = wolfram_str .. "{" else wolfram_str = wolfram_str .. ",{" end
    local str = ""
    for j = 1, #mat do
      str = str .. string.format(" %02d", mat[i][j])
      wolfram_str = wolfram_str .. string.format(j == 1 and "%d" or ",%d", mat[i][j])
    end
    wolfram_str = wolfram_str .. "}"
    print(str)
  end
  wolfram_str = wolfram_str .. "}"
  print(wolfram_str)
end

-- Linear system solver using Cramer's Rule
function DESLOCA.cramer(mat, vec)
  -- Check if matrix is quadratic
	assert(#mat == #mat[1], "Matrix is not square!")
  -- Check if vector has the same size of the matrix
  assert(#mat == #vec, "Vector has not the same size of the matrix!")
	
	local size = #mat
  local main_det = matrix.det(mat)
  local is_indeterminate = false
  if main_det == 0 then is_indeterminate = true end
  
  --[[ DEBUG
  print("\nMat in cramer:")
  print_mat(mat)
  print(string.format("det: %f", main_det))]]
  
  local aux_mats = {}
  local dets = {}
  local result = {}
  for i = 1, size do
    -- Construct the auxiliary matrixes
    aux_mats[i] = matrix.copy(mat)
    for j = 1, size do
      aux_mats[i][j][i] = vec[j]
    end
    
    -- Calculate the auxiliary determinants
    dets[i] = matrix.det(aux_mats[i])
    
    -- Calculate results
    result[i] = dets[i]/main_det
    
    --[[ DEBUG
    print(string.format("\naux_mats[%d]:", i))
    print_mat(aux_mats[i])
    print(string.format("det: %f", matrix.det(aux_mats[i])))]]
  end
  
  --[[ DEBUG
  print("\nResult:")
  local result_str = "{"
  for i = 1, size do
    result_str = result_str .. string.format("%f, ", result[i])
  end
  print(result_str .. "}")]]
  
  return result, is_indeterminate
end

-- Function that calculate the displacements using Cramer's Rule
function DESLOCA.get_displacements()
  DESLOCA.displacements = {}
  
  -- Get displacements
  local result, is_indeterminate = DESLOCA.cramer(DESLOCA.structure_stiff_matrix_simp, DESLOCA.structure_load_vector_simp) -- TODO: Cramer is making some divisions by 0 or really small numbers, seems to be due incompatible structures...
                                                                                                                           -- TODO: Perhaps use a reducted version of the matrix and vector (without the simplified part)
  -- Export indetermination
  DESLOCA.indeterminate_system = is_indeterminate
  
  -- Show indetermination warning message
  indetermination_static_text:Show(DESLOCA.indeterminate_system)
  
  -- Copy from result
  for i = 1, #result do
    DESLOCA.displacements[i] = result[i]
    -- Handle indeterminate system
    if DESLOCA.indeterminate_system then
      DESLOCA.displacements[i] = 0
    end
  end
  
  -- Store displacement info in the bar data
  for i = 1, #DESLOCA.bar do
    -- i
    DESLOCA.bar[i].i.disp_h = DESLOCA.displacements[3*DESLOCA.bar[i].global_node.i - 2]
    DESLOCA.bar[i].i.disp_v = DESLOCA.displacements[3*DESLOCA.bar[i].global_node.i - 1]
    DESLOCA.bar[i].i.disp_m = DESLOCA.displacements[3*DESLOCA.bar[i].global_node.i - 0]
    -- f
    DESLOCA.bar[i].f.disp_h = DESLOCA.displacements[3*DESLOCA.bar[i].global_node.f - 2]
    DESLOCA.bar[i].f.disp_v = DESLOCA.displacements[3*DESLOCA.bar[i].global_node.f - 1]
    DESLOCA.bar[i].f.disp_m = DESLOCA.displacements[3*DESLOCA.bar[i].global_node.f - 0]
    
    -- Store in a single vector that is easier to use (for example in DESLOCA.get_reactions)
    DESLOCA.bar[i].displacements_vector = {DESLOCA.bar[i].i.disp_h, DESLOCA.bar[i].i.disp_v, DESLOCA.bar[i].i.disp_m, DESLOCA.bar[i].f.disp_h, DESLOCA.bar[i].f.disp_v, DESLOCA.bar[i].f.disp_m}
    
    -- DEBUG
    --print(string.format("\nBar %d displacements:\n i: h = %f, v = %f, m = %f \n f: h = %f, v = %f, m = %f", i, DESLOCA.bar[i].i.disp_h, DESLOCA.bar[i].i.disp_v, DESLOCA.bar[i].i.disp_m, DESLOCA.bar[i].f.disp_h, DESLOCA.bar[i].f.disp_v, DESLOCA.bar[i].f.disp_m))
  end
  
end

-- Function that calculate the restraint reactions
function DESLOCA.get_reactions()
  
  -- Init the reactions table
  DESLOCA.reactions = {}
  for j = 1, 3*#DESLOCA.node do
    DESLOCA.reactions[j] = 0
  end
  
  -- Get reactions per bar -- TODO: test if there's a way to get for the whole structure in one go
  for i = 1, #DESLOCA.bar do
  
    DESLOCA.bar[i].reactions_vector = {}
    
    --  Reaction = -V_global + K_global * displacements_vector
    
    for j = 1, 6 do -- vector element
      
      -- Reaction data per bar in a vector form
      DESLOCA.bar[i].reactions_vector[j] = 0
      
      -- equivalent to multiply a matrix with a vector
      for k = 1, 6 do -- sum portion
        
        DESLOCA.bar[i].reactions_vector[j] = DESLOCA.bar[i].reactions_vector[j] + DESLOCA.bar[i].stiff_matrix_global[j][k] * DESLOCA.bar[i].displacements_vector[k]
        
      end
      
      DESLOCA.bar[i].reactions_vector[j] = DESLOCA.bar[i].reactions_vector[j] - DESLOCA.bar[i].load_vector_global[j]
      
      -- Reaction data per bar divided into i and f
      -- i
      if     j == 1 then DESLOCA.bar[i].i.react_h = DESLOCA.bar[i].reactions_vector[j]
      elseif j == 2 then DESLOCA.bar[i].i.react_v = DESLOCA.bar[i].reactions_vector[j]
      elseif j == 3 then DESLOCA.bar[i].i.react_m = DESLOCA.bar[i].reactions_vector[j]
      -- f
      elseif j == 4 then DESLOCA.bar[i].f.react_h = DESLOCA.bar[i].reactions_vector[j]
      elseif j == 5 then DESLOCA.bar[i].f.react_v = DESLOCA.bar[i].reactions_vector[j]
      elseif j == 6 then DESLOCA.bar[i].f.react_m = DESLOCA.bar[i].reactions_vector[j] end
    end
    
    -- Export bar reactions to the structure data
    -- bar i
    DESLOCA.reactions[3*DESLOCA.bar[i].global_node.i - 2] = DESLOCA.reactions[3*DESLOCA.bar[i].global_node.i - 2] + DESLOCA.bar[i].reactions_vector[1]
    DESLOCA.reactions[3*DESLOCA.bar[i].global_node.i - 1] = DESLOCA.reactions[3*DESLOCA.bar[i].global_node.i - 1] + DESLOCA.bar[i].reactions_vector[2]
    DESLOCA.reactions[3*DESLOCA.bar[i].global_node.i - 0] = DESLOCA.reactions[3*DESLOCA.bar[i].global_node.i - 0] + DESLOCA.bar[i].reactions_vector[3]
    -- bar f
    DESLOCA.reactions[3*DESLOCA.bar[i].global_node.f - 2] = DESLOCA.reactions[3*DESLOCA.bar[i].global_node.f - 2] + DESLOCA.bar[i].reactions_vector[4]
    DESLOCA.reactions[3*DESLOCA.bar[i].global_node.f - 1] = DESLOCA.reactions[3*DESLOCA.bar[i].global_node.f - 1] + DESLOCA.bar[i].reactions_vector[5]
    DESLOCA.reactions[3*DESLOCA.bar[i].global_node.f - 0] = DESLOCA.reactions[3*DESLOCA.bar[i].global_node.f - 0] + DESLOCA.bar[i].reactions_vector[6]
    
    -- DEBUG
    --print(string.format("\nBar %d reactions:\n i: h = %f, v = %f, m = %f \n f: h = %f, v = %f, m = %f", i, DESLOCA.bar[i].i.react_h, DESLOCA.bar[i].i.react_v, DESLOCA.bar[i].i.react_m, DESLOCA.bar[i].f.react_h, DESLOCA.bar[i].f.react_v, DESLOCA.bar[i].f.react_m))
  end
  
end

-- Function that calculate the bar forces on its ends in local coordinates
function DESLOCA.get_forces_on_ends()
  
  -- Get forces on ends per bar
  for i = 1, #DESLOCA.bar do
    
    DESLOCA.bar[i].forces_on_ends = {}
    
    local K_R = matrix.mul(DESLOCA.bar[i].stiff_matrix_local, DESLOCA.bar[i].rotation_matrix)
    
    for j = 1, 6 do -- vector element
    
      DESLOCA.bar[i].forces_on_ends[j] = 0
      
      -- equivalent to multiply a matrix with a vector
      for k = 1, 6 do -- sum portion inside vector element
        
        DESLOCA.bar[i].forces_on_ends[j] = DESLOCA.bar[i].forces_on_ends[j] + K_R[j][k] * DESLOCA.bar[i].displacements_vector[k]
        
      end
      
      DESLOCA.bar[i].forces_on_ends[j] = DESLOCA.bar[i].forces_on_ends[j] - DESLOCA.bar[i].load_vector_local[j]
    end
    
    --[[
    local f_global = {
      Rt[1][1]*f_local[1] + Rt[1][2]*f_local[2] + Rt[1][3]*f_local[3] + Rt[1][4]*f_local[4] + Rt[1][5]*f_local[5] + Rt[1][6]*f_local[6],
      Rt[2][1]*f_local[1] + Rt[2][2]*f_local[2] + Rt[2][3]*f_local[3] + Rt[2][4]*f_local[4] + Rt[2][5]*f_local[5] + Rt[2][6]*f_local[6],
      Rt[3][1]*f_local[1] + Rt[3][2]*f_local[2] + Rt[3][3]*f_local[3] + Rt[3][4]*f_local[4] + Rt[3][5]*f_local[5] + Rt[3][6]*f_local[6],
      Rt[4][1]*f_local[1] + Rt[4][2]*f_local[2] + Rt[4][3]*f_local[3] + Rt[4][4]*f_local[4] + Rt[4][5]*f_local[5] + Rt[4][6]*f_local[6],
      Rt[5][1]*f_local[1] + Rt[5][2]*f_local[2] + Rt[5][3]*f_local[3] + Rt[5][4]*f_local[4] + Rt[5][5]*f_local[5] + Rt[5][6]*f_local[6],
      Rt[6][1]*f_local[1] + Rt[6][2]*f_local[2] + Rt[6][3]*f_local[3] + Rt[6][4]*f_local[4] + Rt[6][5]*f_local[5] + Rt[6][6]*f_local[6]
    }]]
  end
end

-- Function that calculate the bar diagram forces
function DESLOCA.get_diagram_forces()
  
  -- Get diagram forces per bar
  for i = 1, #DESLOCA.bar do
    
    -- Get dimensions and angles
    local L, C, S = DESLOCA.dimension(DESLOCA.bar[i])    
    
    -- Init table with default values
    DESLOCA.bar[i].diagram = {N = 0, Qi = 0, Qf = 0, Mi = 0, Mf = 0, Mparabola = false, x_zero = false}
    
    -- Get N
    DESLOCA.bar[i].diagram.N = -DESLOCA.bar[i].forces_on_ends[1] -- could be [4] too, both are x forces and for N it doesn't matter ; "-" because the convention is that compression is negative
    
    -- Get Q
    DESLOCA.bar[i].diagram.Qi = DESLOCA.bar[i].forces_on_ends[2]
    DESLOCA.bar[i].diagram.Qf = -DESLOCA.bar[i].forces_on_ends[5]
    
    -- Get q in local coordinates
    local q_local, q_global
    if DESLOCA.bar[i].dist.is_global then
      q_global = DESLOCA.bar[i].dist.q
      if DESLOCA.bar[i].dist.global_dir == "x" then
        q_local = q_global*S
      elseif DESLOCA.bar[i].dist.global_dir == "y" then
        q_local = q_global*C
      end
    else
      q_local = DESLOCA.bar[i].dist.q
      q_global = q_local*C -- this one is never used   
    end 
    
    -- Get the x position where Q changes signal, if it changes
    if DESLOCA.bar[i].diagram.Qi * DESLOCA.bar[i].diagram.Qf < 0 then -- opposite signals
      DESLOCA.bar[i].diagram.x_zero = math.abs(DESLOCA.bar[i].diagram.Qi/q_local)
    end
    
    -- Get M
    DESLOCA.bar[i].diagram.Mi = -DESLOCA.bar[i].forces_on_ends[3]
    DESLOCA.bar[i].diagram.Mf = DESLOCA.bar[i].forces_on_ends[6]
    
    -- Get Mparabola if Q changes signal
    if DESLOCA.bar[i].diagram.x_zero then
      DESLOCA.bar[i].diagram.Mparabola = DESLOCA.bar[i].diagram.Mi + (DESLOCA.bar[i].diagram.Qi*DESLOCA.bar[i].diagram.x_zero)/2
    end
    
    --[[ DEBUG
    print("\nBarra "..i) -- REMOVE/TEST
    print("N = "..math.round(DESLOCA.bar[i].diagram.N, 1)) -- REMOVE/TEST
    print("Qi = "..math.round(DESLOCA.bar[i].diagram.Qi, 1)) -- REMOVE/TEST
    print("Qf = "..math.round(DESLOCA.bar[i].diagram.Qf, 1)) -- REMOVE/TEST
    print("Mi = "..math.round(DESLOCA.bar[i].diagram.Mi, 1)) -- REMOVE/TEST
    print("Mf = "..math.round(DESLOCA.bar[i].diagram.Mf, 1)) -- REMOVE/TEST
    if DESLOCA.bar[i].diagram.x_zero then
      print("x_zero = "..math.round(DESLOCA.bar[i].diagram.x_zero, 1)) -- REMOVE/TEST
      print("Mparabola = "..math.round(DESLOCA.bar[i].diagram.Mparabola, 1)) -- REMOVE/TEST
    end]]
  end
  
end

-- Icon used, in hex string
DESLOCA.icon_hex = "0000010001004040000001002000284200001600000028000000400000008000000001002000000000000080000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705"..
"921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E"..
"3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FFBBB751FFB"..
"BB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF00000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000BBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FFBBB751FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF958637FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF705921FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF4E3311FF0000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001FFFFFF00000000000FFFFF000000000001FFFF0000000000007FFF0000000000001FFF00000000000007FF00000000000003FF00000000000001FF00000000000000FF000000000000007F000000000000003F000000000000003F000000000000001F000000000000000F000000000000000F0000000000000007000000000000000700000000000000030000000000000003000000000000000300000000000000010000000000000001000000000000000100000000000000010000000000000001000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000100000000000000010000000000000001000000000000000100000000000000030000000000000003000000000000000300000000000000070000000000000007000000000000000F000000000000000F000000000000001F000000000000001F000000000000003F000000000000007F00000000000000FF00000000000001FF00000000000003FF00000000000007FF0000000000001FFF0000000000003FFF000000000001FFFF000000000007FFFF0000000000FFFFFF"


--############################################################################################################################################################
-- WXLUA ENVIROMENT
--############################################################################################################################################################

-- Load wxLua module (needed when not running via wxLua, wxLuaFreeze, or wxLuaEdit)
package.cpath = package.cpath..";./?.dll;./?.so;../lib/?.so;../lib/vc_dll/?.dll;../lib/bcc_dll/?.dll;../lib/mingw_dll/?.dll;"
require("wx")

-- IDs of the controls
local IDs = {
  file_menu_new_project = 2011,

  bar_mode_choice = 0,
  choice_bar = 1,
  
  static_text_L = 2,
  text_ctrl_L = 3,
  static_text_units_L = 4,
  
  static_text_xiyi = 5,
  text_ctrl_xi = 6,
  text_ctrl_yi = 7,
  static_text_units_xiyi = 8,

  static_text_xfxf = 9,
  text_ctrl_xf = 10,
  text_ctrl_yf = 11,
  static_text_units_xfxf = 12,

  static_text_angle = 13,
  spin_ctrl_angle = 14,
  static_text_units_angle = 15,

  static_text_E = 16,
  text_ctrl_E = 17,
  static_text_units_E = 18,

  static_text_A = 19,
  text_ctrl_A = 20,
  static_text_units_A = 21,

  static_text_I = 22,
  text_ctrl_I = 23,
  static_text_units_I = 24,

  check_box_L_fixed = 25,

  create_bar_button = 26,
  delete_bar_button = 27,
  
  dist_load_text_ctrl_q = 28,
  dist_load_option_local = 30,
  dist_load_option_global_x = 1029,
  dist_load_option_global_y = 1030,
  
  nodal_load_static_text_i = 1031,
  nodal_load_text_ctrl_Fx_i = 31,
  nodal_load_text_ctrl_Fy_i = 32,
  nodal_load_text_ctrl_Mz_i = 33,
  nodal_load_option_global_i = 34,
  nodal_load_option_local_i = 35,
  
  nodal_load_static_text_f = 1035,
  nodal_load_text_ctrl_Fx_f = 36,
  nodal_load_text_ctrl_Fy_f = 37,
  nodal_load_text_ctrl_Mz_f = 38,
  nodal_load_option_global_f = 39,
  nodal_load_option_local_f = 40,
  
  mat_estrut_simp_checkbox = 41,
  
  restrain_static_text_i = 1041,
  restrain_option_free_x_i = 42,
  restrain_option_fixed_x_i = 43,
  restrain_option_free_y_i = 44,
  restrain_option_fixed_y_i = 45,
  restrain_option_free_z_i = 46,
  restrain_option_fixed_z_i = 47,
  restrain_static_text_f = 1047,
  restrain_option_free_x_f = 48,
  restrain_option_fixed_x_f = 49,
  restrain_option_free_y_f = 50,
  restrain_option_fixed_y_f = 51,
  restrain_option_free_z_f = 52,
  restrain_option_fixed_z_f = 53,
  
  structure_mode_choice = 1053,
  
  deformed_scale_slider = 1054,
  deformed_scale_static_text = 1055,
  
  diagram_option_N = 2054,
  diagram_option_Q = 2055,
  diagram_option_M = 2056,
  
  grid_results_option_disp = 54,
  grid_results_option_react = 55,
  
  indetermination_static_text = 56,
  
  node_info_static_text_i_1 = 57,
  node_info_static_text_f_1 = 58,
  node_info_static_text_i_2 = 59,
  node_info_static_text_f_2 = 60,
  
  grid_rotation_mat_option_normal = 61,
  grid_rotation_mat_option_transposed = 62,
  
  grid_mat_local = 7001,
  grid_load_vector_local = 7002,
  grid_rotation_mat = 7003,
  grid_mat_global = 7004,
  grid_load_vector_global = 7005,

  grid_f_ep = 7006,
  grid_u = 7007,
  grid_f_local = 7008,

  grid_mat_estrut = 7009,
  grid_load_vector_structure = 7010,
  grid_results = 7011,
  
  grid_load_vector_structure_option1 = 63,
  grid_load_vector_structure_option2 = 64,
  grid_load_vector_structure_option3 = 65,
}

-- Set the default window size
DESLOCA.width, DESLOCA.height = 1280 - 10, 720 - 10 -- -10 to handle the window style correctly

-- Get some system info
local SYSTEM = {}
SYSTEM.screen_width, SYSTEM.screen_height = wx.wxDisplay():GetClientArea():GetWidth(), wx.wxDisplay():GetClientArea():GetHeight()

-- Check if screen is small
if SYSTEM.screen_width < DESLOCA.width or SYSTEM.screen_height < DESLOCA.height then
  SYSTEM.screen_is_small = true
else
  SYSTEM.screen_is_small = false
end

-- Create the window handle
local frame = wx.wxFrame(wx.NULL, wx.wxID_ANY, "DESLOCA", wx.wxDefaultPosition, wx.wxSize(DESLOCA.width, DESLOCA.height),
  wx.wxDEFAULT_FRAME_STYLE - (wx.wxRESIZE_BORDER + wx.wxMAXIMIZE_BOX)) -- a style combination that makes the window not resizable
frame:Centre(wx.wxBOTH) -- to center the window horizontally and vertically

-- Create the menu handles
local fileMenu = wx.wxMenu("", wx.wxMENU_TEAROFF)
fileMenu:Append(IDs.file_menu_new_project, "Novo projeto", "Iniciar um projeto em branco")
fileMenu:Append(wx.wxID_EXIT, "Sair", "Fechar o DESLOCA")

local optionsMenu = wx.wxMenu()

local showMenu = wx.wxMenu()

local helpMenu = wx.wxMenu()
helpMenu:Append(wx.wxID_ABOUT, "Sobre", "Sobre o DESLOCA")

-- Create a menu bar and append the menus
local menuBar = wx.wxMenuBar()
menuBar:Append(fileMenu, "&Arquivo")
--menuBar:Append(optionsMenu, "&Opções") -- TODO: implement later
--menuBar:Append(showMenu, "&Exibir") -- TODO: implement later
menuBar:Append(helpMenu, "&Ajuda")

frame:SetMenuBar(menuBar)

-- Create a status bar to display messages
status_bar = frame:CreateStatusBar(3) -- 3 fields
local status_bar_field_widths = {250, 300, 1000}
status_bar:SetStatusWidths(status_bar_field_widths)

-- Statically indeterminate structure warning
status_bar:SetForegroundColour(wx.wxRED)
indetermination_static_text = wx.wxStaticText(status_bar, wx.wxID_ANY, "Estrutura estaticamente indeterminada! Verificar vinculações ou conexões entre as barras!",
                                              wx.wxPoint(status_bar_field_widths[1] + status_bar_field_widths[2] + 8, 2), wx.wxDefaultSize)
indetermination_static_text:SetFont(wx.wxFont(12, wx.wxFONTFAMILY_DEFAULT, wx.wxFONTSTYLE_NORMAL, wx.wxNORMAL))

-- Create a single child window, wxWidgets will set the size to fill frame
panel = wx.wxPanel(frame, wx.wxID_ANY, wx.wxPoint(20, 0), wx.wxSize(200,500))

-- Rename some functions
local fmt = string.format

-- Rounding function, num is number and dec is decimal places to round (from http://lua-users.org/wiki/SimpleRound by Igor Skoric)
function math.round(num, dec)
  local mult = 10^(dec or 0)
  if num >= 0 then
    return math.floor(num * mult + 0.5) / mult
  else
    return math.ceil(num * mult - 0.5) / mult
  end
end

--############################################################################################################################################################
-- PAINT
--############################################################################################################################################################

-- Painting functions and constants
local Paint = {}

-- Initialization of some constants
Paint.bar_sel = 1

Paint.canvas_x, Paint.canvas_y = 4, 345
Paint.canvas_width, Paint.canvas_height = 344, 280 --380, 280
Paint.canvas_end_x, Paint.canvas_end_y = Paint.canvas_x + Paint.canvas_width, Paint.canvas_y + Paint.canvas_height
Paint.canvas_scale = 0.20 -- 1 means each 1*100 pixels is 1 meter of the structure

Paint.bar_canvas_x, Paint.bar_canvas_y = 4, 25
Paint.bar_canvas_width, Paint.bar_canvas_height = 258, 190
Paint.bar_canvas_end_x, Paint.bar_canvas_end_y = Paint.bar_canvas_x + Paint.bar_canvas_width, Paint.bar_canvas_y + Paint.bar_canvas_height

Paint.forces_canvas_x, Paint.forces_canvas_y = Paint.bar_canvas_end_x + 230, Paint.bar_canvas_y + 160 -- +230+398 to snap it in the right border
Paint.forces_canvas_width, Paint.forces_canvas_height = 370, 124
Paint.forces_canvas_end_x, Paint.forces_canvas_end_y = Paint.forces_canvas_x + Paint.forces_canvas_width, Paint.forces_canvas_y + Paint.forces_canvas_height

-- Paint event handler for the frame that's called by wxEVT_PAINT
function Paint.paint(event)

  -- Create the DC handle, in order to draw stuff
  Paint.dc = wx.wxPaintDC(panel) -- wxBufferedPaintDC or wxAutoBufferedPaintDC to avoid flickering, but background is black
  
  -- Change default font to a smaller one
  Paint.dc:SetFont(wx.wxFont(8, wx.wxFONTFAMILY_DEFAULT, wx.wxFONTSTYLE_NORMAL, wx.wxFONTWEIGHT_NORMAL)) -- TODO: figure out a way to restore the default font
  
  -- SUB-FUNCTIONS --
  
  -- Draw an arrow given (x1, y1), (x2, y2) and a size
  local function draw_arrow(x1, y1, x2, y2, head)
    
    local angle = math.atan((y2-y1)/(x2-x1)) -- in radians -- TODO: perhaps use the method from DESLOCA.dimension()
    local head_size = head or 10
    
    local angle1, angle2 = angle + math.pi/4, angle - math.pi/4 --0.785398163398, angle - 0.785398163398 -- 45° in radians
    local delta_x1, delta_y1 = head_size*math.cos(angle1), head_size*math.sin(angle1)
    local delta_x2, delta_y2 = head_size*math.cos(angle2), head_size*math.sin(angle2)
    local head1_x1, head1_y1 = x2, y2 
    local head1_x2, head1_y2 
    local head2_x1, head2_y1 = x2, y2
    local head2_x2, head2_y2
    
    if x1 < x2 then -- 1st and 4th quadrant
      head1_x2, head1_y2 = head1_x1 - delta_x1, head1_y1 - delta_y1
      head2_x2, head2_y2 = head2_x1 - delta_x2, head2_y1 - delta_y2
    elseif x1 == x2 then -- vertical arrow
      head1_x2, head1_y2 = head1_x1 - delta_x1, head1_y1 - delta_y1
      head2_x2, head2_y2 = head2_x1 - delta_x2, head2_y1 - delta_y2
    else
      head1_x2, head1_y2 = head1_x1 + delta_x1, head1_y1 + delta_y1
      head2_x2, head2_y2 = head2_x1 + delta_x2, head2_y1 + delta_y2
    end
    
    -- Draw
    Paint.dc:DrawLine(x1, y1, x2, y2)
    Paint.dc:DrawLine(head1_x1, head1_y1, head1_x2, head1_y2)
    Paint.dc:DrawLine(head2_x1, head2_y1, head2_x2, head2_y2)
  end
  
  -- Restraint drawning  
  local function draw_restraint(type, x_pos, y_pos, free_z)
    Paint.dc:SetPen(wx.wxPen(wx.wxLIGHT_GREY, 2, 1))
    
    local node_x, node_y = x_pos, y_pos
    local diff = free_z and 0 or 8 -- to avoid overlap of triangles and squares
    
    if type == "x" then -- horizontal triangle
      Paint.dc:DrawLine(node_x - diff, node_y, node_x - 14 - diff, node_y - 8)
      Paint.dc:DrawLine(node_x - diff, node_y, node_x - 14 - diff, node_y + 8)
      Paint.dc:DrawLine(node_x - 14 - diff, node_y - 8, node_x - 14 - diff, node_y + 8)
    elseif type == "y" then -- vertical triangle
      Paint.dc:DrawLine(node_x, node_y + diff, node_x - 8, node_y + 14 + diff)
      Paint.dc:DrawLine(node_x, node_y + diff, node_x + 8, node_y + 14 + diff)
      Paint.dc:DrawLine(node_x - 8, node_y + 14 + diff, node_x + 8, node_y + 14 + diff)
    elseif type == "z" then -- square
      Paint.dc:DrawLine(node_x - 8, node_y - 8, node_x + 8, node_y - 8)
      Paint.dc:DrawLine(node_x - 8, node_y - 8, node_x - 8, node_y + 8)
      Paint.dc:DrawLine(node_x + 8, node_y + 8, node_x + 8, node_y - 8)
      Paint.dc:DrawLine(node_x + 8, node_y + 8, node_x - 8, node_y + 8)
    end
  end
  
  -- Distributed load drawing
  local function draw_dist_load(bar, xi, yi, xf, yf)
    Paint.dc:SetPen(wx.wxPen(wx.wxColour("LIME GREEN"), 1, 1))
    Paint.dc:SetTextForeground(wx.wxColour("LIME GREEN"))
    
    -- Import bar properties
    local L, C, S = DESLOCA.dimension(bar)
    local angle_rad, angle_deg = Get_angle(S, C)
    local q, is_global = bar.dist.q, bar.dist.is_global
    
    -- Don't draw if q = 0
    if q == 0 then return end
    
    -- Signal
    local signal = q < 0 and -1 or 1
    
    -- Arrow constants
    local arrow_n = math.floor(math.round(L, 1)/100) -- round to avoid, for example, L being 99,99999 instead of 100 sometimes due trigonometry roundings
    if arrow_n > 10 then arrow_n = 10 -- to avoid too many arrows in the bar canvas
    elseif arrow_n < 4 then arrow_n = 4 end -- to avoid too few arrows in both canvas
    local arrow_dx, arrow_dy = (xf - xi)/(arrow_n - 1), (yf - yi)/(arrow_n - 1)
    
    -- q string
    local q_str = fmt("%.1f kN/cm", math.round(signal*q, 1))
    local q_str_w, q_str_h = Paint.dc:GetTextExtent(q_str)
    q_str_w, q_str_h = q_str_w - 20, q_str_h/2 -- because GetTextExtent is not correct for DrawRotatedText
    local q_str_x = (xi + xf)/2 - C*q_str_w/2 + signal*(20 + 2*q_str_h)*math.cos(angle_rad - math.pi/2) -- TODO: correct for positive q
    local q_str_y = (yi + yf)/2 + S*q_str_w/2 - signal*(20 + 2*q_str_h)*math.sin(angle_rad - math.pi/2) -- TODO: correct for positive q
    if angle_deg >= 90 and angle_deg < 270 then
      angle_deg = (angle_deg + 180)%360 -- to avoid upside-down text
      q_str_x = (xi + xf)/2 + C*q_str_w/2 + signal*(20 + 0*q_str_h)*math.cos(angle_rad - math.pi/2) -- to correct its position accordingly -- TODO: correct for positive q
      q_str_y = (yi + yf)/2 - S*q_str_w/2 - signal*(20 + 0*q_str_h)*math.sin(angle_rad - math.pi/2) -- to correct its position accordingly -- TODO: correct for positive q
    end
    
    -- Drawing
    if is_global then -- global q
      
      -- in x
      if bar.dist.global_dir == "x" then
        Paint.dc:DrawLine(xi - signal*20, yi, xf - signal*20, yf)
        for i = 0, arrow_n - 1 do
          draw_arrow(xi + i*arrow_dx - signal*20, yi + i*arrow_dy, xi + i*arrow_dx, yi + i*arrow_dy, 6)
        end
        Paint.dc:DrawRotatedText(q_str, (xi + xf)/2 - signal*40, (yi + yf)/2, angle_deg)
      -- in y
      elseif bar.dist.global_dir == "y" then
        Paint.dc:DrawLine(xi, yi + signal*20, xf, yf + signal*20)
        for i = 0, arrow_n - 1 do
          draw_arrow(xi + i*arrow_dx, yi + i*arrow_dy + signal*20, xi + i*arrow_dx, yi + i*arrow_dy, 6)
        end
        Paint.dc:DrawRotatedText(q_str, (xi + xf)/2, (yi + yf)/2 + signal*40, angle_deg)
      end
      
    else -- local q
    
      Paint.dc:DrawLine(xi + signal*20*math.cos(angle_rad - math.pi/2), yi - signal*20*math.sin(angle_rad - math.pi/2), xf + signal*20*math.cos(angle_rad - math.pi/2), yf - signal*20*math.sin(angle_rad - math.pi/2))
      for i = 0, arrow_n - 1 do
        draw_arrow(xi + signal*20*math.cos(angle_rad - math.pi/2) + i*arrow_dx, yi - signal*20*math.sin(angle_rad - math.pi/2) + i*arrow_dy, xi + i*arrow_dx, yi + i*arrow_dy, 6)
      end
      Paint.dc:DrawRotatedText(q_str, q_str_x, q_str_y, angle_deg)
      
    end
  end
  
  -- Nodal load drawings
  local function draw_nodal_load(xi, yi, xf, yf, canvas, bar, node)
    Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x00, 0xA4, 0xA4), 2, 1)) -- petrol green
    Paint.dc:SetTextForeground(wx.wxColour(0x00, 0xA4, 0xA4)) -- petrol green
    Paint.dc:SetBrush(wx.wxBrush(wx.wxNullColour, wx.wxTRANSPARENT))
    
    if canvas == "bar" then -- nodal load drawing in bar canvas, uses loads related to the bar
    
      -- Import bar properties
      local L, C, S = DESLOCA.dimension(bar)
      local angle_rad, angle_deg = Get_angle(S, C)
      local Fx_i, Fy_i, Mz_i, is_global_i = bar.i.Fx, bar.i.Fy, bar.i.Mz, bar.i.is_global
      local Fx_f, Fy_f, Mz_f, is_global_f = bar.f.Fx, bar.f.Fy, bar.f.Mz, bar.f.is_global
      
      -- Get F values to be printed
      local F_i = math.sqrt(Fx_i^2 + Fy_i^2)
      local F_f = math.sqrt(Fx_f^2 + Fy_f^2)
      
      -- Mz signals
      local Mz_i_signal = Mz_i < 0 and -1 or 1
      local Mz_f_signal = Mz_f < 0 and -1 or 1
      
      -- i
      if F_i ~= 0 then
        
        -- Get the load arrow angles based on Fx and Fy
        local C_i, S_i = Fx_i/F_i, Fy_i/F_i
        local angle_rad_i, angle_deg_i = Get_angle(S_i, C_i)
        
        -- F_i string
        local F_i_str = fmt("%.1f kN", math.round(F_i, 1))
        local F_i_str_w, F_i_str_h = Paint.dc:GetTextExtent(F_i_str)
        F_i_str_w, F_i_str_h = F_i_str_w - 14, F_i_str_h/2 -- because GetTextExtent is not correct for DrawRotatedText
        
        if is_global_i then
        
          draw_arrow(xi - 40*C_i, yi + 40*S_i, xi - 10*C_i, yi + 10*S_i, 8)
          Paint.dc:DrawRotatedText(fmt("%.1f kN", math.round(F_i, 1)), xi - (40 + 2*F_i_str_h)*C_i + (F_i_str_w/2)*math.cos(angle_rad_i - math.pi/2),
                                                        yi + (40 + 2*F_i_str_h)*S_i - 2*F_i_str_h*math.sin(angle_rad_i - math.pi/2), angle_deg_i + 90)
        else
          
          draw_arrow(xi - 40*math.cos(angle_rad + angle_rad_i), yi + 40*math.sin(angle_rad + angle_rad_i),
                     xi - 10*math.cos(angle_rad + angle_rad_i), yi + 10*math.sin(angle_rad + angle_rad_i), 8)
          Paint.dc:DrawRotatedText(fmt("%.1f kN", math.round(F_i, 1)),xi - (40 + 2*F_i_str_h)*math.cos(angle_rad + angle_rad_i) + (F_i_str_w/2)*math.cos(angle_rad + angle_rad_i - math.pi/2),
                                                       yi + (40 + 2*F_i_str_h)*math.sin(angle_rad + angle_rad_i) - 2*F_i_str_h*math.sin(angle_rad + angle_rad_i - math.pi/2), angle_deg + angle_deg_i + 90)
        end
      end
      if Mz_i ~= 0 then
        -- Arrow constants
        local arrow_r = 20
        local arrow_x, arrow_y = xi - 0.86603*arrow_r*Mz_i_signal, yi - 0.5*arrow_r -- cos(30°) = 0.86603, sin(30°) = 0.5
        
        -- Draw the load
        Paint.dc:DrawEllipticArc(xi - arrow_r, yi - arrow_r, 2*arrow_r, 2*arrow_r, 30, 150)
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x - 8*0.25882*Mz_i_signal, arrow_y - 8*0.96593) -- sin(15°) = 0.25882, cos(15°) = 0.96593
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x + 8*0.96593*Mz_i_signal, arrow_y - 8*0.25882)  
        Paint.dc:DrawText(fmt("%.1f kNcm", math.round(Mz_i*Mz_i_signal, 1)), xi - Paint.dc:GetTextExtent(fmt("%.1f kNcm", Mz_i*Mz_i_signal))/2, yi - arrow_r - 16)
      end
      
      -- f
      if F_f ~= 0 then
      
        -- Get the load arrow angles based on Fx and Fy
        local C_f, S_f = Fx_f/F_f, Fy_f/F_f
        local angle_rad_f, angle_deg_f = Get_angle(S_f, C_f)
        
        -- F_f string
        local F_f_str = fmt("%.1f kN", math.round(F_f, 1))
        local F_f_str_w, F_f_str_h = Paint.dc:GetTextExtent(F_f_str)
        F_f_str_w, F_f_str_h = F_f_str_w - 14, F_f_str_h/2 -- because GetTextExtent is not correct for DrawRotatedText
        
        if is_global_f then
        
          draw_arrow(xf - 40*C_f, yf + 40*S_f, xf - 10*C_f, yf + 10*S_f, 8)
          Paint.dc:DrawRotatedText(fmt("%.1f kN", math.round(F_f, 1)), xf - (40 + 2*F_f_str_h)*C_f + (F_f_str_w/2)*math.cos(angle_rad_f - math.pi/2),
                                                        yf + (40 + 2*F_f_str_h)*S_f - 2*F_f_str_h*math.sin(angle_rad_f - math.pi/2), angle_deg_f + 90)
        else
          
          draw_arrow(xf - 40*math.cos(angle_rad + angle_rad_f), yf + 40*math.sin(angle_rad + angle_rad_f),
                     xf - 10*math.cos(angle_rad + angle_rad_f), yf + 10*math.sin(angle_rad + angle_rad_f), 8)
          Paint.dc:DrawRotatedText(fmt("%.1f kN", math.round(F_f, 1)),xf - (40 + 2*F_f_str_h)*math.cos(angle_rad + angle_rad_f) + (F_f_str_w/2)*math.cos(angle_rad + angle_rad_f - math.pi/2),
                                                       yf + (40 + 2*F_f_str_h)*math.sin(angle_rad + angle_rad_f) - 2*F_f_str_h*math.sin(angle_rad + angle_rad_f - math.pi/2), angle_deg + angle_deg_f + 90)
        end
      end
      if Mz_f ~= 0 then
        -- Arrow constants
        local arrow_r = 20
        local arrow_x, arrow_y = xf - 0.86603*arrow_r*Mz_f_signal, yf - 0.5*arrow_r -- cos(30°) = 0.86603, sin(30°) = 0.5
        
        -- Draw the load
        Paint.dc:DrawEllipticArc(xf - arrow_r, yf - arrow_r, 2*arrow_r, 2*arrow_r, 30, 150)
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x - 8*0.25882*Mz_f_signal, arrow_y - 8*0.96593) -- sin(15°) = 0.25882, cos(15°) = 0.96593
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x + 8*0.96593*Mz_f_signal, arrow_y - 8*0.25882)  
        Paint.dc:DrawText(fmt("%.1f kNcm", math.round(Mz_f*Mz_f_signal, 1)), xf - Paint.dc:GetTextExtent(fmt("%.1f kNcm", Mz_f*Mz_f_signal))/2, yf - arrow_r - 16)
      end
      
    elseif canvas == "structure" then -- nodal load drawing in structure canvas, uses loads related to the node
      
      local x, y = xi, yi
      
      -- Import node loads
      local Fx = node.Fx
      local Fy = node.Fy
      local Mz = node.Mz
      
      -- Get F value to be printed
      local F = math.sqrt(Fx^2 + Fy^2)
      
      -- Print F if not null
      if F ~= 0 then
        -- Get the load arrow angles based on Fx and Fy
        local C, S = Fx/F, Fy/F
        local angle_rad, angle_deg = Get_angle(S, C)
        
        -- F string
        local F_str = fmt("%.1f kN", math.round(F, 1))
        local F_str_w, F_str_h = Paint.dc:GetTextExtent(F_str)
        F_str_w, F_str_h = F_str_w - 14, F_str_h/2 -- because GetTextExtent is not correct for DrawRotatedText
        
        -- Draw the load
        draw_arrow(x - 40*C, y + 40*S, x - 10*C, y + 10*S, 8)
        Paint.dc:DrawRotatedText(fmt("%.1f kN", math.round(F, 1)), x - (40 + 2*F_str_h)*C + (F_str_w/2)*math.cos(angle_rad - math.pi/2),
                                                                   y + (40 + 2*F_str_h)*S - 2*F_str_h*math.sin(angle_rad - math.pi/2), angle_deg + 90)
      end
      
      -- Print Mz if not null
      if Mz ~= 0 then
        -- Get the load signal, to properly draw the arrow
        local signal = Mz < 0 and -1 or 1
        
        -- Arrow constants
        local arrow_r = 20
        local arrow_x, arrow_y = x - 0.86603*arrow_r*signal, y - 0.5*arrow_r -- cos(30°) = 0.86603, sin(30°) = 0.5
        
        -- Draw the load
        Paint.dc:DrawEllipticArc(x - arrow_r, y - arrow_r, 2*arrow_r, 2*arrow_r, 30, 150)
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x - 8*0.25882*signal, arrow_y - 8*0.96593) -- sin(15°) = 0.25882, cos(15°) = 0.96593
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x + 8*0.96593*signal, arrow_y - 8*0.25882)  
        Paint.dc:DrawText(fmt("%.1f kNcm", math.round(Mz*signal, 1)), x - Paint.dc:GetTextExtent(fmt("%.1f kNcm", Mz*signal))/2, y - arrow_r - 16)
      end
      
    end
  end
  
  -- Nodal reactions drawing
  local function draw_nodal_reactions(canvas, x, y, bar, direction)
    Paint.dc:SetTextForeground(wx.wxColour(0xFF, 0x80, 0x00)) -- orange
    Paint.dc:SetPen(wx.wxPen(wx.wxColour(0xFF, 0x80, 0x00), 2, 1)) -- orange
    Paint.dc:SetBrush(wx.wxBrush(wx.wxNullColour, wx.wxTRANSPARENT))
    Paint.dc:SetFont(wx.wxFont(9, wx.wxFONTFAMILY_DEFAULT, wx.wxFONTSTYLE_NORMAL, wx.wxFONTWEIGHT_NORMAL))
    
    if canvas == "structure" then
      
      for i = 1, #DESLOCA.node do
      
        local draw_x, draw_y = x + Paint.canvas_scale*DESLOCA.node[i].x, y - Paint.canvas_scale*DESLOCA.node[i].y
        local reaction_value, signal
        
        if not DESLOCA.node[i].free_x then
          reaction_value = math.round(DESLOCA.reactions[i*3 - 2], 2)
          signal = reaction_value < 0 and -1 or 1
          
          draw_arrow(draw_x - 40*signal, draw_y, draw_x - 10*signal, draw_y, 8)
          
          local x_pos
          if signal == 1 then
            x_pos = draw_x - 10 - Paint.dc:GetTextExtent(fmt("%.2f kN", math.round(reaction_value*signal, 2)))
          else
            x_pos = draw_x + 20
          end
          Paint.dc:DrawText(fmt("%.2f kN", math.round(reaction_value*signal, 2)), x_pos, draw_y + 4)
        end
        if not DESLOCA.node[i].free_y then
          reaction_value = math.round(DESLOCA.reactions[i*3 - 1], 2)
          signal = reaction_value < 0 and -1 or 1
        
          draw_arrow(draw_x, draw_y + 25 + 15*signal, draw_x, draw_y + 25 - 15*signal, 8)
          Paint.dc:DrawText(fmt("%.2f kN", math.round(reaction_value*signal, 2)), draw_x + 2, draw_y + 40)          
        end
        if not DESLOCA.node[i].free_z then
          reaction_value = math.round(DESLOCA.reactions[i*3 - 0], 2)
          signal = reaction_value < 0 and -1 or 1
        
          local arrow_r = 25
          local arrow_x, arrow_y = draw_x - 0.86603*arrow_r*signal, draw_y - 0.5*arrow_r -- cos(30°) = 0.86603, sin(30°) = 0.5
          Paint.dc:DrawEllipticArc(draw_x - arrow_r, draw_y - arrow_r, 2*arrow_r, 2*arrow_r, 30, 150)
          Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x - 8*0.25882*signal, arrow_y - 8*0.96593) -- sin(15°) = 0.25882, cos(15°) = 0.96593
          Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x + 8*0.96593*signal, arrow_y - 8*0.25882)  
          Paint.dc:DrawText(fmt("%.2f kNcm", math.round(reaction_value*signal, 2)), draw_x - Paint.dc:GetTextExtent(fmt("%.2f kNcm", reaction_value*signal))/2, draw_y - arrow_r - 16)
        end
      end
      
    elseif canvas == "forces on ends" then
      
      local reaction_value = bar.forces_on_ends[direction]
      local signal = reaction_value < 0 and -1 or 1
      local arrow_r = 25 -- for z only
      local arrow_x, arrow_y = x - 0.86603*arrow_r*signal, y - 0.5*arrow_r -- for z only; cos(30°) = 0.86603, sin(30°) = 0.5
      
      if     direction == 1 then -- i x
        draw_arrow(x - 25 - 15*signal, y, x - 25 + 15*signal, y, 8)
        Paint.dc:DrawText(fmt("%.2f kN", math.round(reaction_value*signal, 2)), x - 10 - Paint.dc:GetTextExtent(fmt("%.2f kN", reaction_value*signal)), y + 6)
      elseif direction == 2 then -- i y
        draw_arrow(x, y + 25 + 15*signal, x, y + 25 - 15*signal, 8)
        Paint.dc:DrawText(fmt("%.2f kN", math.round(reaction_value*signal, 2)), x - Paint.dc:GetTextExtent(fmt("%.2f kN", reaction_value*signal))/2, y + 41)  
      elseif direction == 3 then -- i z
        Paint.dc:DrawEllipticArc(x - arrow_r, y - arrow_r, 2*arrow_r, 2*arrow_r, 30, 150)
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x - 8*0.25882*signal, arrow_y - 8*0.96593) -- sin(15°) = 0.25882, cos(15°) = 0.96593
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x + 8*0.96593*signal, arrow_y - 8*0.25882)
        Paint.dc:DrawText(fmt("%.2f kNcm", math.round(reaction_value*signal, 2)), x - Paint.dc:GetTextExtent(fmt("%.2f kNcm", reaction_value*signal))/2, y - arrow_r - 16)
      elseif direction == 4 then -- f x
        draw_arrow(x + 25 - 15*signal, y, x + 25 + 15*signal, y, 8)
        Paint.dc:DrawText(fmt("%.2f kN", math.round(reaction_value*signal, 2)), x + 10, y + 6)
      elseif direction == 5 then -- f y
        draw_arrow(x, y + 25 + 15*signal, x, y + 25 - 15*signal, 8)
        Paint.dc:DrawText(fmt("%.2f kN", math.round(reaction_value*signal, 2)), x - Paint.dc:GetTextExtent(fmt("%.2f kN", reaction_value*signal))/2, y + 41) 
      elseif direction == 6 then -- f z
        Paint.dc:DrawEllipticArc(x - arrow_r, y - arrow_r, 2*arrow_r, 2*arrow_r, 30, 150)
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x - 8*0.25882*signal, arrow_y - 8*0.96593) -- sin(15°) = 0.25882, cos(15°) = 0.96593
        Paint.dc:DrawLine(arrow_x, arrow_y, arrow_x + 8*0.96593*signal, arrow_y - 8*0.25882)
        Paint.dc:DrawText(fmt("%.2f kNcm", math.round(reaction_value*signal, 2)), x - Paint.dc:GetTextExtent(fmt("%.2f kNcm", reaction_value*signal))/2, y - arrow_r - 16)
      end
      
    end
  end
  
  --- STRUCTURE ---
  
  -- Measure the structure
  local smaller_x, smaller_y, bigger_x, bigger_y = 1000000, 1000000, -1000000, -1000000 -- arbitrary, could've used math.huge
  for i = 1, #DESLOCA.bar do
    if DESLOCA.bar[i].xi < smaller_x then smaller_x = DESLOCA.bar[i].xi end
    if DESLOCA.bar[i].xf < smaller_x then smaller_x = DESLOCA.bar[i].xf end
    if DESLOCA.bar[i].yi < smaller_y then smaller_y = DESLOCA.bar[i].yi end
    if DESLOCA.bar[i].yf < smaller_y then smaller_y = DESLOCA.bar[i].yf end
    
    if DESLOCA.bar[i].xi > bigger_x then bigger_x = DESLOCA.bar[i].xi end
    if DESLOCA.bar[i].xf > bigger_x then bigger_x = DESLOCA.bar[i].xf end
    if DESLOCA.bar[i].yi > bigger_y then bigger_y = DESLOCA.bar[i].yi end
    if DESLOCA.bar[i].yf > bigger_y then bigger_y = DESLOCA.bar[i].yf end
  end
  Paint.center_x = math.floor(Paint.canvas_scale*(bigger_x + smaller_x)/2)
  Paint.center_y = math.floor(Paint.canvas_scale*(bigger_y + smaller_y)/2)
  
  -- Draw the structure grid
  Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x69, 0x76, 0x99), 1, 1)) -- ice
  Paint.dc:SetBrush(wx.wxBrush(wx.wxColour(0x21, 0x28, 0x30), wx.wxSOLID)) -- peacot blue (AutoCAD background)
  Paint.dc:DrawRectangle(Paint.canvas_x, Paint.canvas_y, Paint.canvas_width, Paint.canvas_height)
  
  local origin_x, origin_y = Paint.canvas_x + Paint.canvas_width/2, Paint.canvas_y + Paint.canvas_height/2

  local _, x_dec = math.modf((origin_x - Paint.center_x - Paint.canvas_x)/(Paint.canvas_scale*100)) -- to get how "misaligned" the (origin_x-Paint.center_x) is related with Paint.canvas_x
  local _, y_dec = math.modf((origin_y + Paint.center_y - Paint.canvas_y)/(Paint.canvas_scale*100)) -- to get how "misaligned" the (origin_y-Paint.center_y) is related with Paint.canvas_y
  local grid_start_x = Paint.canvas_x + x_dec*(Paint.canvas_scale*100) -- top left grid dot, "aligned" with (origin_x-Paint.center_x)
  local grid_start_y = Paint.canvas_y + y_dec*(Paint.canvas_scale*100) -- top left grid dot, "aligned" with (origin_y-Paint.center_y)
  
  for i = grid_start_x, grid_start_x + Paint.canvas_width, Paint.canvas_scale*100 do
    for j = grid_start_y, grid_start_y + Paint.canvas_height, Paint.canvas_scale*100 do
      Paint.dc:DrawPoint(i, j)
    end
  end
  
  -- Draw restraints
  for i = 1, #DESLOCA.node do
    local draw_x, draw_y = origin_x + Paint.canvas_scale*DESLOCA.node[i].x - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.node[i].y + Paint.center_y

    if DESLOCA.structure_mode ~= 1 and DESLOCA.structure_mode ~= 3 then -- not "Reações" nor "Diagramas de esforços internos"
      if not DESLOCA.node[i].free_x then draw_restraint("x", draw_x, draw_y, DESLOCA.node[i].free_z) end
      if not DESLOCA.node[i].free_y then draw_restraint("y", draw_x, draw_y, DESLOCA.node[i].free_z) end
      if not DESLOCA.node[i].free_z then draw_restraint("z", draw_x, draw_y) end
    end
  end
  
  -- Draw the structure
  for i = 1, #DESLOCA.bar do
    -- DEBUG
    --print("Bar " .. i) -- REMOVE/TEST
  
    local draw_xi, draw_yi = origin_x + Paint.canvas_scale*DESLOCA.bar[i].xi - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.bar[i].yi + Paint.center_y
    local draw_xf, draw_yf = origin_x + Paint.canvas_scale*DESLOCA.bar[i].xf - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.bar[i].yf + Paint.center_y
    
    -- bars
    if i == Paint.bar_sel then -- highlight selected bar
      Paint.dc:SetPen(wx.wxPen(wx.wxColour("ORANGE"), 5, 1))
    else
      Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x00, 0x40, 0xFF), 3, 1)) -- blue
      --Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x7F, 0x59, 0x3F), 3, 1)) -- light brown 7F593F
      --Paint.dc:SetPen(wx.wxPen(wx.wxColour(0xFF, 0xE5, 0xB4), 3, 1)) -- peach FFE5B4
      --Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x7F, 0xA1, 0xFF), 3, 1)) -- light blue 7FA1FF
    end
    Paint.dc:DrawLine(draw_xi, draw_yi, draw_xf, draw_yf)
    
    if DESLOCA.structure_mode == 0 then -- "Estrutura com cargas"
      -- bar ids
      Paint.dc:SetBrush(wx.wxWHITE_BRUSH)
      if i == Paint.bar_sel then -- highlight selected bar
        Paint.dc:SetPen(wx.wxPen(wx.wxColour("ORANGE"), 1, 1))
        Paint.dc:SetTextForeground(wx.wxColour("ORANGE"))
      else
        Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x00, 0x40, 0xFF), 1, 1)) -- blue
        Paint.dc:SetTextForeground(wx.wxColour(0x00, 0x40, 0xFF)) -- blue
      end
      
      local id_str = fmt("%d", i)
      local id_str_w, id_str_h = Paint.dc:GetTextExtent(id_str)
      Paint.dc:DrawRectangle(origin_x + Paint.canvas_scale*(DESLOCA.bar[i].xi + DESLOCA.bar[i].xf)/2 - Paint.center_x - id_str_w/2 - 3,
                             origin_y - Paint.canvas_scale*(DESLOCA.bar[i].yi + DESLOCA.bar[i].yf)/2 + Paint.center_y - id_str_h + 5,
                             id_str_w + 6, id_str_h + 3)
      Paint.dc:DrawText(id_str, origin_x + Paint.canvas_scale*(DESLOCA.bar[i].xi + DESLOCA.bar[i].xf)/2 - Paint.center_x - id_str_w/2,
                                origin_y - Paint.canvas_scale*(DESLOCA.bar[i].yi + DESLOCA.bar[i].yf)/2 + Paint.center_y - id_str_h + 6
      )
      
      -- distributed load
      draw_dist_load(DESLOCA.bar[i], draw_xi, draw_yi, draw_xf, draw_yf)
    end
  end
  
  -- Draw nodal loads
  for i = 1, #DESLOCA.node do
    local draw_x, draw_y = origin_x + Paint.canvas_scale*DESLOCA.node[i].x - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.node[i].y + Paint.center_y
    
    if DESLOCA.structure_mode == 0 then -- "Estrutura com cargas"
      draw_nodal_load(draw_x, draw_y, nil, nil, "structure", nil, DESLOCA.node[i])
    end
  end
  
  -- Draw the nodes, with labels
  Paint.dc:SetTextForeground(wx.wxColour(0xFF, 0xD4, 0x00)) -- yellow
  Paint.dc:SetPen(wx.wxPen(wx.wxColour(0xFF, 0xD4, 0x00), 3, 1)) -- yellow
  for i = 1, #DESLOCA.node do
    Paint.dc:DrawCircle(origin_x + Paint.canvas_scale*DESLOCA.node[i].x - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.node[i].y + Paint.center_y, 2)
    Paint.dc:DrawText(fmt("%d", i), origin_x + Paint.canvas_scale*DESLOCA.node[i].x - Paint.center_x + 6, origin_y - Paint.canvas_scale*DESLOCA.node[i].y + Paint.center_y + 6)
  end

  -- Draw the reactions
  if DESLOCA.structure_mode == 1 then -- "Reações"
    draw_nodal_reactions("structure", origin_x - Paint.center_x, origin_y + Paint.center_y)
  end
  
  -- Draw the deformed configuration
  if DESLOCA.structure_mode == 2 then -- "Configuração deformada"
    Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x00, 0xFF, 0xFF), 1, 1)) -- cyan
    
    -- Scan thru the bars to get the greatest deformation
    local max_def = 0
    for j = 1, #DESLOCA.bar do
      -- nodal displacements
      if math.abs(DESLOCA.bar[j].displacements_vector[1]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[1]) end -- i dx
      if math.abs(DESLOCA.bar[j].displacements_vector[2]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[2]) end -- i dy
      if math.abs(DESLOCA.bar[j].displacements_vector[4]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[4]) end -- f dx
      if math.abs(DESLOCA.bar[j].displacements_vector[5]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[5]) end -- f dy
      
      local L = DESLOCA.dimension(DESLOCA.bar[j])
      -- peak deformation due distributed load
      local dist_load_peak_def = DESLOCA.bar[j].dist.q * L^4 / (384 * DESLOCA.bar[j].E * DESLOCA.bar[j].I) -- peak deformation q*L^4/(384*E*I), from w(x) = q*x^4/(24*E*I) with x = L/2
      if math.abs(dist_load_peak_def) > max_def then max_def = math.abs(dist_load_peak_def) end
      
      --[[ DEBUG
      print("\nBar "..j) -- REMOVE/TEST
      print("i dx: "..DESLOCA.bar[j].displacements_vector[1]) -- REMOVE/TEST
      print("i dy: "..DESLOCA.bar[j].displacements_vector[2]) -- REMOVE/TEST
      print("f dx: "..DESLOCA.bar[j].displacements_vector[4]) -- REMOVE/TEST
      print("f dy: "..DESLOCA.bar[j].displacements_vector[5]) -- REMOVE/TEST
      print("dist_load_peak_def: "..dist_load_peak_def) -- REMOVE/TEST]]
    end
    if max_def == 0 then max_def = 1 end -- to avoid diving by 0 later
    --DEBUG
    --print("\nmax_def: "..max_def) -- REMOVE/TEST
    
    -- Constants used
    local max_deformed_scale = 300/max_def
    local deformed_scale = (max_deformed_scale/4)*deformed_scale_slider:GetValue()
    local rot_size = 75 -- used in "i/f rotation" to force the spline simulating the property that bars mantain their relative angles in a rotation (i.e they doesn't have hinges)
    
    -- Updates the scale number static text accordingly
    deformed_scale_static_text:SetLabel(fmt("Escala\n%.1f", math.round(deformed_scale, 1)))
    
    -- Handle indetermination
    if DESLOCA.indeterminate_system then deformed_scale = 0 end
    
    -- DEBUG
    --print("max_deformed_scale: "..max_deformed_scale) -- REMOVE/TEST
    --print("deformed_scale: "..deformed_scale) -- REMOVE/TEST
    
    -- Draw the deformed configuration for each bar
    for j = 1, #DESLOCA.bar do
      
      local draw_xi, draw_yi = origin_x + Paint.canvas_scale*DESLOCA.bar[j].xi - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.bar[j].yi + Paint.center_y
      local draw_xf, draw_yf = origin_x + Paint.canvas_scale*DESLOCA.bar[j].xf - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.bar[j].yf + Paint.center_y
      
      local L, C, S = DESLOCA.dimension(DESLOCA.bar[j])
      local angle_rad = Get_angle(S, C)
      
      local deformed_angle_i = angle_rad + DESLOCA.bar[j].displacements_vector[3] * deformed_scale/4 -- dividing by 4 to avoid making it too deformed
      local deformed_angle_f = angle_rad + math.pi + DESLOCA.bar[j].displacements_vector[6] * deformed_scale/4 -- dividing by 4 to avoid making it too deformed
      
      -- Get q in local coordinates
      local q_local, q_global
      if DESLOCA.bar[j].dist.is_global then
        q_global = DESLOCA.bar[j].dist.q
        if DESLOCA.bar[j].dist.global_dir == "x" then
          q_local = q_global*S
        elseif DESLOCA.bar[j].dist.global_dir == "y" then
          q_local = q_global*C
        end
      else
        q_local = DESLOCA.bar[j].dist.q
        q_global = q_local*C -- this one is never used   
      end 
      local dist_load_peak_def = q_local * L^4 / (384 * DESLOCA.bar[j].E * DESLOCA.bar[j].I) -- peak deformation q*L^4/(384*E*I), from w(x) = q*x^4/(24*E*I) with x = L/2
      
      local deformed_coords = {[1] = {}, [2] = {}, [3] = {}, [4] = {}, [5] = {}}
      -- i
      deformed_coords[1][1] = draw_xi + Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[1] * deformed_scale
      deformed_coords[1][2] = draw_yi - Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[2] * deformed_scale
      -- i rotation
      deformed_coords[2][1] = draw_xi + Paint.canvas_scale*rot_size*math.cos(deformed_angle_i) + Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[1] * deformed_scale
      deformed_coords[2][2] = draw_yi - Paint.canvas_scale*rot_size*math.sin(deformed_angle_i) - Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[2] * deformed_scale
      -- f rotation
      deformed_coords[4][1] = draw_xf + Paint.canvas_scale*rot_size*math.cos(deformed_angle_f) + Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[4] * deformed_scale
      deformed_coords[4][2] = draw_yf - Paint.canvas_scale*rot_size*math.sin(deformed_angle_f) - Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[5] * deformed_scale
      -- f
      deformed_coords[5][1] = draw_xf + Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[4] * deformed_scale
      deformed_coords[5][2] = draw_yf - Paint.canvas_scale*DESLOCA.bar[j].displacements_vector[5] * deformed_scale
      -- bar middle (deformation due distributed load)
      deformed_coords[3][1] = (deformed_coords[1][1] + deformed_coords[5][1])/2 - Paint.canvas_scale*dist_load_peak_def*math.sin(angle_rad) * deformed_scale
      deformed_coords[3][2] = (deformed_coords[1][2] + deformed_coords[5][2])/2 - Paint.canvas_scale*dist_load_peak_def*math.cos(angle_rad) * deformed_scale
      
      Paint.dc:DrawSpline(deformed_coords)
      
    end
  end
  
  -- Draw the force diagrams
  if DESLOCA.structure_mode == 3 and not DESLOCA.indeterminate_system then -- "Diagramas de esforços internos"
    Paint.dc:SetTextForeground(wx.wxColour(0x00, 0xFF, 0xFF)) -- cyan
    
    -- Scan thru the bars to get the greatest N, Q and M
    local max_N, max_Q, max_M = 0, 0, 0
    for i = 1, #DESLOCA.bar do
      -- N
      if math.abs(DESLOCA.bar[i].diagram.N) > max_N then max_N = math.abs(DESLOCA.bar[i].diagram.N) end
      -- Q
      if math.abs(DESLOCA.bar[i].diagram.Qi) > max_Q then max_Q = math.abs(DESLOCA.bar[i].diagram.Qi) end
      if math.abs(DESLOCA.bar[i].diagram.Qf) > max_Q then max_Q = math.abs(DESLOCA.bar[i].diagram.Qf) end
      -- M
      if math.abs(DESLOCA.bar[i].diagram.Mi) > max_M then max_M = math.abs(DESLOCA.bar[i].diagram.Mi) end
      if math.abs(DESLOCA.bar[i].diagram.Mf) > max_M then max_M = math.abs(DESLOCA.bar[i].diagram.Mf) end
      if DESLOCA.bar[i].diagram.x_zero then
        if math.abs(DESLOCA.bar[i].diagram.Mparabola) > max_M then max_M = math.abs(DESLOCA.bar[i].diagram.Mparabola) end
      end
    end    
    
    -- Draw the diagram for the selected force type, one bar at a time
    for i = 1, #DESLOCA.bar do
      Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x00, 0xFF, 0xFF), 1, 1)) -- cyan
      
      -- Get correct drawing coordinates for nodes i and f
      local draw_xi, draw_yi = origin_x + Paint.canvas_scale*DESLOCA.bar[i].xi - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.bar[i].yi + Paint.center_y
      local draw_xf, draw_yf = origin_x + Paint.canvas_scale*DESLOCA.bar[i].xf - Paint.center_x, origin_y - Paint.canvas_scale*DESLOCA.bar[i].yf + Paint.center_y
      
      -- Get dimensions
      local L, C, S = DESLOCA.dimension(DESLOCA.bar[i])
      local angle_rad, angle_deg = Get_angle(S, C)
    
      -- Get q in local coordinates
      local q_local, q_global
      if DESLOCA.bar[i].dist.is_global then
        q_global = DESLOCA.bar[i].dist.q
        if DESLOCA.bar[i].dist.global_dir == "x" then
          q_local = q_global*S
        elseif DESLOCA.bar[i].dist.global_dir == "y" then
          q_local = q_global*C
        end
      else
        q_local = DESLOCA.bar[i].dist.q
        q_global = q_local*C -- this one is never used   
      end 
      
      local diagram_scale
      
      if diagram_option_N:GetValue() == true then -- N
        diagram_scale = 300/max_N
        
        Paint.dc:DrawLine(draw_xi, draw_yi, draw_xi - DESLOCA.bar[i].diagram.N * S *(diagram_scale*Paint.canvas_scale), draw_yi - DESLOCA.bar[i].diagram.N * C *(diagram_scale*Paint.canvas_scale))
        Paint.dc:DrawLine(draw_xf, draw_yf, draw_xf - DESLOCA.bar[i].diagram.N * S *(diagram_scale*Paint.canvas_scale), draw_yf - DESLOCA.bar[i].diagram.N * C *(diagram_scale*Paint.canvas_scale))
        Paint.dc:DrawLine(draw_xi - DESLOCA.bar[i].diagram.N * S *(diagram_scale*Paint.canvas_scale), draw_yi - DESLOCA.bar[i].diagram.N * C *(diagram_scale*Paint.canvas_scale),
                          draw_xf - DESLOCA.bar[i].diagram.N * S *(diagram_scale*Paint.canvas_scale), draw_yf - DESLOCA.bar[i].diagram.N * C *(diagram_scale*Paint.canvas_scale)
        )
        
        -- Signal
        local N_signal = DESLOCA.bar[i].diagram.N < 0 and -1 or 1
        
        -- N string
        local N_str = fmt("%.1f", math.round(DESLOCA.bar[i].diagram.N, 1))
        local N_str_w, N_str_h = Paint.dc:GetTextExtent(N_str)
        
        Paint.dc:DrawRotatedText(N_str,
          (draw_xi + draw_xf)/2 - (DESLOCA.bar[i].diagram.N) * S *(diagram_scale*Paint.canvas_scale) - N_str_w/2 * C - (N_signal == 1 and (N_str_h * S) or 0),
          (draw_yi + draw_yf)/2 - (DESLOCA.bar[i].diagram.N) * C *(diagram_scale*Paint.canvas_scale) + N_str_w/2 * S - (N_signal == 1 and (N_str_h * C) or 0),
        angle_deg)
        
      elseif diagram_option_Q:GetValue() == true then -- Q
        diagram_scale = 300/max_Q
        
        Paint.dc:DrawLine(draw_xi, draw_yi, draw_xi - DESLOCA.bar[i].diagram.Qi * S *(diagram_scale*Paint.canvas_scale), draw_yi - DESLOCA.bar[i].diagram.Qi * C *(diagram_scale*Paint.canvas_scale))
        Paint.dc:DrawLine(draw_xf, draw_yf, draw_xf - DESLOCA.bar[i].diagram.Qf * S *(diagram_scale*Paint.canvas_scale), draw_yf - DESLOCA.bar[i].diagram.Qf * C *(diagram_scale*Paint.canvas_scale))
        Paint.dc:DrawLine(draw_xi - DESLOCA.bar[i].diagram.Qi * S *(diagram_scale*Paint.canvas_scale), draw_yi - DESLOCA.bar[i].diagram.Qi * C *(diagram_scale*Paint.canvas_scale),
                          draw_xf - DESLOCA.bar[i].diagram.Qf * S *(diagram_scale*Paint.canvas_scale), draw_yf - DESLOCA.bar[i].diagram.Qf * C *(diagram_scale*Paint.canvas_scale)
        )
        
        -- Signals
        local Qi_signal = DESLOCA.bar[i].diagram.Qi < 0 and -1 or 1
        local Qf_signal = DESLOCA.bar[i].diagram.Qf < 0 and -1 or 1
        
        -- Qi string
        local Qi_str = fmt("%.1f", math.round(DESLOCA.bar[i].diagram.Qi, 1))
        local Qi_str_w, Qi_str_h = Paint.dc:GetTextExtent(Qi_str)
        --Qi_str_h = Qi_str_h - 5 -- because GetTextExtent is not correct for DrawRotatedText
        
        -- Qf string
        local Qf_str = fmt("%.1f", math.round(DESLOCA.bar[i].diagram.Qf, 1))
        local Qf_str_w, Qf_str_h = Paint.dc:GetTextExtent(Qf_str)
        --Qf_str_h = Qf_str_h - 5 -- because GetTextExtent is not correct for DrawRotatedText
        
        Paint.dc:DrawRotatedText(Qi_str,
          draw_xi - DESLOCA.bar[i].diagram.Qi * S *(diagram_scale*Paint.canvas_scale) - (Qi_signal == 1 and Qi_str_h * S or 0),
          draw_yi - DESLOCA.bar[i].diagram.Qi * C *(diagram_scale*Paint.canvas_scale) - (Qi_signal == 1 and Qi_str_h * C or 0),
        angle_deg)
        Paint.dc:DrawRotatedText(Qf_str,
          draw_xf - DESLOCA.bar[i].diagram.Qf * S *(diagram_scale*Paint.canvas_scale) - Qf_str_w * C - (Qf_signal == 1 and Qf_str_h * S or 0),
          draw_yf - DESLOCA.bar[i].diagram.Qf * C *(diagram_scale*Paint.canvas_scale) + Qf_str_w * S - (Qf_signal == 1 and Qf_str_h * C or 0),
        angle_deg)
        
      elseif diagram_option_M:GetValue() == true then -- M
        diagram_scale = 300/max_M
        
        Paint.dc:DrawLine(draw_xi, draw_yi, draw_xi + DESLOCA.bar[i].diagram.Mi * S *(diagram_scale*Paint.canvas_scale), draw_yi + DESLOCA.bar[i].diagram.Mi * C *(diagram_scale*Paint.canvas_scale))
        Paint.dc:DrawLine(draw_xf, draw_yf, draw_xf + DESLOCA.bar[i].diagram.Mf * S *(diagram_scale*Paint.canvas_scale), draw_yf + DESLOCA.bar[i].diagram.Mf * C *(diagram_scale*Paint.canvas_scale))
        
        -- Signals
        local Mi_signal = DESLOCA.bar[i].diagram.Mi < 0 and -1 or 1
        local Mf_signal = DESLOCA.bar[i].diagram.Mf < 0 and -1 or 1
        
        -- Mi string
        local Mi_str = fmt("%.1f", math.round(DESLOCA.bar[i].diagram.Mi, 1))
        local Mi_str_w, Mi_str_h = Paint.dc:GetTextExtent(Mi_str)
        --Mi_str_h = Mi_str_h - 5 -- because GetTextExtent is not correct for DrawRotatedText
        
        -- Mf string
        local Mf_str = fmt("%.1f", math.round(DESLOCA.bar[i].diagram.Mf, 1))
        local Mf_str_w, Mf_str_h = Paint.dc:GetTextExtent(Mf_str)
        --Mf_str_h = Mf_str_h - 5 -- because GetTextExtent is not correct for DrawRotatedText
        
        Paint.dc:DrawRotatedText(Mi_str,
          draw_xi + DESLOCA.bar[i].diagram.Mi * S *(diagram_scale*Paint.canvas_scale) - (Mi_signal == -1 and Mi_str_h * S or 0),
          draw_yi + DESLOCA.bar[i].diagram.Mi * C *(diagram_scale*Paint.canvas_scale) - (Mi_signal == -1 and Mi_str_h * C or 0),
        angle_deg)
        Paint.dc:DrawRotatedText(Mf_str,
          draw_xf + DESLOCA.bar[i].diagram.Mf * S *(diagram_scale*Paint.canvas_scale) - Mf_str_w * C - (Mf_signal == -1 and Mf_str_h * S or 0),
          draw_yf + DESLOCA.bar[i].diagram.Mf * C *(diagram_scale*Paint.canvas_scale) + Mf_str_w * S - (Mf_signal == -1 and Mf_str_h * C or 0),
        angle_deg)
        
        -- Draw M curve
        local step = math.round(L, 1)/10 -- round to avoid, for example, L being 99,99999 instead of 100 sometimes due trigonometry roundings
        local M_prev, x_pos_prev, y_pos_prev = DESLOCA.bar[i].diagram.Mi, draw_xi + DESLOCA.bar[i].diagram.Mi * S *(diagram_scale*Paint.canvas_scale), draw_yi + DESLOCA.bar[i].diagram.Mi * C *(diagram_scale*Paint.canvas_scale)
        local Q, M, x_pos, y_pos
        for x = step, L+0.001, step do -- starting at step to not start at zero ; L+0.001 to make sure the last step is drawn, otherwise sometimes it's not due roundings
          
          Q = DESLOCA.bar[i].diagram.Qi + q_local * x -- gerenal Q equation
          
          M = DESLOCA.bar[i].diagram.Mi + ((DESLOCA.bar[i].diagram.Qi + Q) * x)/2 -- general M equation
          
          x_pos = draw_xi + x * C * Paint.canvas_scale + M * S *(diagram_scale*Paint.canvas_scale)
          y_pos = draw_yi - x * S * Paint.canvas_scale + M * C *(diagram_scale*Paint.canvas_scale)
          
          Paint.dc:DrawLine( x_pos_prev, y_pos_prev, x_pos, y_pos)
          
          M_prev, x_pos_prev, y_pos_prev = M, x_pos, y_pos
        end
        
        -- Parabola max, if it exists, indicated with a dashed line
        if DESLOCA.bar[i].diagram.x_zero then
          local x_zero_x_pos = draw_xi + Paint.canvas_scale * DESLOCA.bar[i].diagram.x_zero * C
          local x_zero_y_pos = draw_yi - Paint.canvas_scale * DESLOCA.bar[i].diagram.x_zero * S
          
          -- Mparabola signal
          local Mp_signal = DESLOCA.bar[i].diagram.Mparabola < 0 and -1 or 1
        
          -- Mparabola string
          local Mp_str = fmt("%.1f", math.round(DESLOCA.bar[i].diagram.Mparabola, 1))
          local Mp_str_w, Mp_str_h = Paint.dc:GetTextExtent(Mp_str)
          --Mp_str_h = Mp_str_h - 5 -- because GetTextExtent is not correct for DrawRotatedText
          
          Paint.dc:SetPen(wx.wxPen(wx.wxColour(0x00, 0xFF, 0xFF), 1, wx.wxDOT)) -- cyan ; wxDOT is a dashed line, for some reason
          Paint.dc:DrawLine(x_zero_x_pos, x_zero_y_pos, x_zero_x_pos + DESLOCA.bar[i].diagram.Mparabola * S *(diagram_scale*Paint.canvas_scale), x_zero_y_pos + DESLOCA.bar[i].diagram.Mparabola * C *(diagram_scale*Paint.canvas_scale))
        
          Paint.dc:DrawRotatedText(Mp_str,
            x_zero_x_pos + DESLOCA.bar[i].diagram.Mparabola * S *(diagram_scale*Paint.canvas_scale) - Mp_str_w/2 * C - (Mp_signal == -1 and Mp_str_h * S or 0),
            x_zero_y_pos + DESLOCA.bar[i].diagram.Mparabola * C *(diagram_scale*Paint.canvas_scale) + Mp_str_w/2 * S - (Mp_signal == -1 and Mp_str_h * C or 0),
          angle_deg)
        end
        
      end
    
    end
    
    --[[
    -- Scan thru the bars to get the greatest deformation
    local max_def = 0
    for j = 1, #DESLOCA.bar do
      -- nodal displacements
      if math.abs(DESLOCA.bar[j].displacements_vector[1]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[1]) end -- i dx
      if math.abs(DESLOCA.bar[j].displacements_vector[2]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[2]) end -- i dy
      if math.abs(DESLOCA.bar[j].displacements_vector[4]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[4]) end -- f dx
      if math.abs(DESLOCA.bar[j].displacements_vector[5]) > max_def then max_def = math.abs(DESLOCA.bar[j].displacements_vector[5]) end -- f dy
      
      local L = DESLOCA.dimension(DESLOCA.bar[j])
      -- peak deformation due distributed load
      local dist_load_peak_def = DESLOCA.bar[j].dist.q * L^4 / (384 * DESLOCA.bar[j].E * DESLOCA.bar[j].I) -- peak deformation q*L^4/(384*E*I), from w(x) = q*x^4/(24*E*I) with x = L/2
      if math.abs(dist_load_peak_def) > max_def then max_def = math.abs(dist_load_peak_def) end
      
      -- DEBUG
      print("\nBar "..j) -- REMOVE/TEST
      print("i dx: "..DESLOCA.bar[j].displacements_vector[1]) -- REMOVE/TEST
      print("i dy: "..DESLOCA.bar[j].displacements_vector[2]) -- REMOVE/TEST
      print("f dx: "..DESLOCA.bar[j].displacements_vector[4]) -- REMOVE/TEST
      print("f dy: "..DESLOCA.bar[j].displacements_vector[5]) -- REMOVE/TEST
      print("dist_load_peak_def: "..dist_load_peak_def) -- REMOVE/TEST
    end
    if max_def == 0 then max_def = 1 end -- to avoid diving by 0 later]]
  end
  
  --- BACKGROUND (to "avoid" drawings leaking out from the canvases) --- TODO: remove if structure/bar size limit is created
  
  -- Redraw canvas border to prioritize it
  Paint.dc:SetPen(wx.wxBLACK_PEN)
  Paint.dc:DrawLine(Paint.canvas_x, Paint.canvas_y, Paint.canvas_x, Paint.canvas_end_y-1)
  Paint.dc:DrawLine(Paint.canvas_x, Paint.canvas_y, Paint.canvas_end_x-1, Paint.canvas_y)
  Paint.dc:DrawLine(Paint.canvas_x, Paint.canvas_end_y-1, Paint.canvas_end_x-1, Paint.canvas_end_y-1)
  Paint.dc:DrawLine(Paint.canvas_end_x-1, Paint.canvas_y, Paint.canvas_end_x-1, Paint.canvas_end_y)
  
  -- Draw background with the panel colour
  Paint.dc:SetPen(wx.wxPen(panel:GetBackgroundColour(), 1, 1))
  Paint.dc:SetBrush(wx.wxBrush(panel:GetBackgroundColour(), wx.wxSOLID))
  Paint.dc:DrawRectangle(0, 0, DESLOCA.width, Paint.canvas_y)
  Paint.dc:DrawRectangle(0, Paint.canvas_end_y, DESLOCA.width, DESLOCA.height - Paint.canvas_end_y)
  Paint.dc:DrawRectangle(0, Paint.canvas_y, Paint.canvas_x, Paint.canvas_height)
  Paint.dc:DrawRectangle(Paint.canvas_end_x, Paint.canvas_y, DESLOCA.width - Paint.canvas_end_x, Paint.canvas_height)
  
  --- SELECTED BAR ---
  
  -- DEBUG
  --print("Selected bar " .. Paint.bar_sel) -- REMOVE/TEST
  Paint.dc:SetFont(wx.wxFont(8, wx.wxFONTFAMILY_DEFAULT, wx.wxFONTSTYLE_NORMAL, wx.wxFONTWEIGHT_NORMAL))
  
  -- Draw the bar canvas
  Paint.dc:SetPen(wx.wxBLACK_PEN)
  Paint.dc:SetBrush(wx.wxBrush(wx.wxColour(0x21, 0x28, 0x30), wx.wxSOLID)) -- peacot blue (AutoCAD background)
  Paint.dc:DrawRectangle(Paint.bar_canvas_x, Paint.bar_canvas_y, Paint.bar_canvas_width, Paint.bar_canvas_height)
  
  -- Drawing constants
  local L, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
  L = Paint.bar_canvas_height - 100 -- to fit the bar and its info inside the canvas
  local angle_rad, angle_deg = Get_angle(S, C)
  local dx, dy = C*L/2, S*L/2
  Paint.bar_center_x, Paint.bar_center_y = Paint.bar_canvas_x + Paint.bar_canvas_width/2, Paint.bar_canvas_y + Paint.bar_canvas_height/2
  
  -- Restraints
  if DESLOCA.bar[Paint.bar_sel].i.free_x == false then draw_restraint("x", Paint.bar_center_x - dx, Paint.bar_center_y + dy, DESLOCA.bar[Paint.bar_sel].i.free_z) end
  if DESLOCA.bar[Paint.bar_sel].i.free_y == false then draw_restraint("y", Paint.bar_center_x - dx, Paint.bar_center_y + dy, DESLOCA.bar[Paint.bar_sel].i.free_z) end
  if DESLOCA.bar[Paint.bar_sel].i.free_z == false then draw_restraint("z", Paint.bar_center_x - dx, Paint.bar_center_y + dy) end
  if DESLOCA.bar[Paint.bar_sel].f.free_x == false then draw_restraint("x", Paint.bar_center_x + dx, Paint.bar_center_y - dy, DESLOCA.bar[Paint.bar_sel].f.free_z) end
  if DESLOCA.bar[Paint.bar_sel].f.free_y == false then draw_restraint("y", Paint.bar_center_x + dx, Paint.bar_center_y - dy, DESLOCA.bar[Paint.bar_sel].f.free_z) end
  if DESLOCA.bar[Paint.bar_sel].f.free_z == false then draw_restraint("z", Paint.bar_center_x + dx, Paint.bar_center_y - dy) end
  
  -- Draw the bar
  Paint.dc:SetPen(wx.wxPen(wx.wxColour("ORANGE"), 5, 1))
  Paint.dc:DrawLine(Paint.bar_center_x - dx, Paint.bar_center_y + dy, Paint.bar_center_x + dx, Paint.bar_center_y - dy)
  
  -- Draw the nodes
  Paint.dc:SetPen(wx.wxPen(wx.wxColour(0xFF, 0xD4, 0x00), 3, 1)) -- yellow
  Paint.dc:DrawCircle(Paint.bar_center_x - dx, Paint.bar_center_y + dy, 2)
  Paint.dc:DrawCircle(Paint.bar_center_x + dx, Paint.bar_center_y - dy, 2)
  
  -- Draw the node names
  Paint.dc:SetPen(wx.wxPen(wx.wxWHITE, 1, wx.wxDOT)) -- wxDOT is a dashed line, for some reason
  Paint.dc:SetTextForeground(wx.wxWHITE)
  -- i
  Paint.dc:DrawLine(Paint.bar_center_x - dx, Paint.bar_center_y + dy, Paint.bar_center_x - dx + 40*math.cos(angle_rad - math.pi/4), Paint.bar_center_y + dy - 40*math.sin(angle_rad - math.pi/4))
  Paint.dc:DrawCircle(Paint.bar_center_x - dx + 40*math.cos(angle_rad - math.pi/4), Paint.bar_center_y + dy - 40*math.sin(angle_rad - math.pi/4), 8)
  Paint.dc:DrawText("i", Paint.bar_center_x - dx + 40*math.cos(angle_rad - math.pi/4), Paint.bar_center_y + dy - 40*math.sin(angle_rad - math.pi/4) - 7)
  Paint.dc:DrawText(fmt("(%d)", DESLOCA.bar[Paint.bar_sel].global_node.i), Paint.bar_center_x - dx + 40*math.cos(angle_rad - math.pi/4) - 2 + 10, Paint.bar_center_y + dy - 40*math.sin(angle_rad - math.pi/4) - 8) -- TODO: make this better
  -- f
  Paint.dc:DrawLine(Paint.bar_center_x + dx, Paint.bar_center_y - dy, Paint.bar_center_x + dx + 40*math.cos(angle_rad - 3*math.pi/4), Paint.bar_center_y - dy - 40*math.sin(angle_rad - 3*math.pi/4))
  Paint.dc:DrawCircle(Paint.bar_center_x + dx + 40*math.cos(angle_rad - 3*math.pi/4), Paint.bar_center_y - dy - 40*math.sin(angle_rad - 3*math.pi/4), 8)
  Paint.dc:DrawText("f", Paint.bar_center_x + dx + 40*math.cos(angle_rad - 3*math.pi/4) - 2, Paint.bar_center_y - dy - 40*math.sin(angle_rad - 3*math.pi/4) - 7)
  Paint.dc:DrawText(fmt("(%d)", DESLOCA.bar[Paint.bar_sel].global_node.f), Paint.bar_center_x + dx + 40*math.cos(angle_rad - 3*math.pi/4) - 2 + 10, Paint.bar_center_y - dy - 40*math.sin(angle_rad - 3*math.pi/4) - 8) -- TODO: make this better
  
  -- Distributed load
  draw_dist_load(DESLOCA.bar[Paint.bar_sel], Paint.bar_center_x - dx, Paint.bar_center_y + dy, Paint.bar_center_x + dx, Paint.bar_center_y - dy)
  
  -- Nodal loads
  draw_nodal_load(Paint.bar_center_x - dx, Paint.bar_center_y + dy, Paint.bar_center_x + dx, Paint.bar_center_y - dy, "bar", DESLOCA.bar[Paint.bar_sel])
  
  --- SELECTED BAR FORCES IN ENDS ---
  
  if bar_mode_choice:GetCurrentSelection() == 1 then -- "Forças nas extremidades em coordenadas locais"
    
    local draw_xi, draw_yi = Paint.forces_canvas_x + Paint.forces_canvas_width/2 - 80, Paint.forces_canvas_y + Paint.forces_canvas_height/2 - 8
    local draw_xf, draw_yf = Paint.forces_canvas_x + Paint.forces_canvas_width/2 + 80, Paint.forces_canvas_y + Paint.forces_canvas_height/2 - 8
    
    -- Draw the canvas
    Paint.dc:SetPen(wx.wxBLACK_PEN)
    Paint.dc:SetBrush(wx.wxBrush(wx.wxColour(0x21, 0x28, 0x30), wx.wxSOLID)) -- peacot blue (AutoCAD background)
    Paint.dc:DrawRectangle(Paint.forces_canvas_x, Paint.forces_canvas_y, Paint.forces_canvas_width, Paint.forces_canvas_height)
    
    -- Draw the bar
    Paint.dc:SetPen(wx.wxPen(wx.wxColour("ORANGE"), 5, 1))
    Paint.dc:DrawLine(draw_xi, draw_yi, draw_xf, draw_yf)
    
    -- Local axis
    Paint.dc:SetPen(wx.wxWHITE_PEN)
    Paint.dc:SetTextForeground(wx.wxWHITE)
    local axis_size = 20
    Paint.dc:DrawText("local", Paint.forces_canvas_x + 10, Paint.forces_canvas_end_y - 20)
    draw_arrow(Paint.forces_canvas_x + 10, Paint.forces_canvas_end_y - 10 - 15, Paint.forces_canvas_x + 10 + axis_size, Paint.forces_canvas_end_y - 10 - 15, 8) -- x
    Paint.dc:DrawText("x", Paint.forces_canvas_x + 10 + axis_size + 5, Paint.forces_canvas_end_y - 10 - 10 - 13) -- x
    draw_arrow(Paint.forces_canvas_x + 10, Paint.forces_canvas_end_y - 10 - 15, Paint.forces_canvas_x + 10, Paint.forces_canvas_end_y - 10 - axis_size - 15, 8) -- y
    Paint.dc:DrawText("y", Paint.forces_canvas_x + 10 - 2, Paint.forces_canvas_end_y - 10 - axis_size - 32) -- y
    
    -- Draw the forces
    -- i
    draw_nodal_reactions("forces on ends", draw_xi, draw_yi, DESLOCA.bar[Paint.bar_sel], 1)
    draw_nodal_reactions("forces on ends", draw_xi, draw_yi, DESLOCA.bar[Paint.bar_sel], 2)
    draw_nodal_reactions("forces on ends", draw_xi, draw_yi, DESLOCA.bar[Paint.bar_sel], 3)
    -- f
    draw_nodal_reactions("forces on ends", draw_xf, draw_yf, DESLOCA.bar[Paint.bar_sel], 4)
    draw_nodal_reactions("forces on ends", draw_xf, draw_yf, DESLOCA.bar[Paint.bar_sel], 5)
    draw_nodal_reactions("forces on ends", draw_xf, draw_yf, DESLOCA.bar[Paint.bar_sel], 6)
    
    -- Draw the operation symbols in the equation
    Paint.dc:SetPen(wx.wxPen(wx.wxBLACK, 2, 1))
    local sym_x, sym_y = Paint.bar_canvas_x + 547, Paint.bar_canvas_y + 85
    -- +
    Paint.dc:DrawLine(sym_x, sym_y, sym_x + 10, sym_y)
    Paint.dc:DrawLine(sym_x + 6, sym_y - 5, sym_x + 6, sym_y + 4) ; sym_x = sym_x + 326 + 8
    -- .
    Paint.dc:DrawCircle(sym_x + 2, sym_y, 2) ; sym_x = sym_x + 206 + 2
    -- .
    Paint.dc:DrawCircle(sym_x + 2, sym_y, 2) ; sym_x = sym_x + 69 + 2
    -- =
    Paint.dc:DrawLine(sym_x, sym_y - 2, sym_x + 10, sym_y - 2)
    Paint.dc:DrawLine(sym_x, sym_y + 2, sym_x + 10, sym_y + 2)
  end
  
  --- MISC ---
  
  Paint.dc:SetPen(wx.wxWHITE_PEN)
  Paint.dc:SetTextForeground(wx.wxWHITE)
  Paint.dc:SetFont(wx.wxFont(8, wx.wxFONTFAMILY_DEFAULT, wx.wxFONTSTYLE_NORMAL, wx.wxFONTWEIGHT_NORMAL))
  
  -- Global axis in structure canvas
  local axis_size = 30
  draw_arrow(Paint.canvas_x + 10, Paint.canvas_end_y - 10, Paint.canvas_x + 10 + axis_size, Paint.canvas_end_y - 10) -- x
  Paint.dc:DrawText("x", Paint.canvas_x + 10 + axis_size + 5, Paint.canvas_end_y - 10 - 8)
  draw_arrow(Paint.canvas_x + 10, Paint.canvas_end_y - 10, Paint.canvas_x + 10, Paint.canvas_end_y - 10 - axis_size) -- y
  Paint.dc:DrawText("y", Paint.canvas_x + 10 - 2, Paint.canvas_end_y - 10 - axis_size - 18)
  
  -- Global axis in bar canvas
  axis_size = 20
  Paint.dc:DrawText("global", Paint.bar_canvas_x + 10 - 3, Paint.bar_canvas_end_y - 19)
  draw_arrow(Paint.bar_canvas_x + 10, Paint.bar_canvas_end_y - 10 - 15, Paint.bar_canvas_x + 10 + axis_size, Paint.bar_canvas_end_y - 10 - 15, 8) -- x
  Paint.dc:DrawText("x", Paint.bar_canvas_x + 10 + axis_size + 5, Paint.bar_canvas_end_y - 10 - 10 - 13) -- x
  draw_arrow(Paint.bar_canvas_x + 10, Paint.bar_canvas_end_y - 10 - 15, Paint.bar_canvas_x + 10, Paint.bar_canvas_end_y - 10 - axis_size - 15, 8) -- y
  Paint.dc:DrawText("y", Paint.bar_canvas_x + 10 - 2, Paint.bar_canvas_end_y - 10 - axis_size - 33) -- y
  -- Local axis in bar canvas
  local spin_r = 15*math.sqrt(2)
  local axis_x, axis_y = Paint.bar_canvas_end_x - 25 - spin_r*math.cos(angle_rad+math.pi/4), Paint.bar_canvas_end_y - 45 + spin_r*math.sin(angle_rad+math.pi/4)
  Paint.dc:DrawText("local", Paint.bar_canvas_end_x - 30, Paint.bar_canvas_end_y - 19)
  draw_arrow(axis_x, axis_y, axis_x + axis_size*C, axis_y - axis_size*S, 8) -- x
  Paint.dc:DrawText("x", axis_x + (axis_size + 6)*C - 2, axis_y - (axis_size + 6)*S - 8) -- x
  draw_arrow(axis_x, axis_y, axis_x - axis_size*S, axis_y - axis_size*C, 8) -- y
  Paint.dc:DrawText("y", axis_x - (axis_size + 8)*S - 2, axis_y - (axis_size + 8)*C - 8) -- y
  
  -- Grid size info
  Paint.dc:DrawText("grid: 100cm", Paint.canvas_end_x - Paint.dc:GetTextExtent("grid: 100cm") - 5, Paint.canvas_end_y - 20)
  
  -- Coordinates system vertical labels for the bar matrixes
  if bar_mode_choice:GetCurrentSelection() == 0 then -- "Matrizes de rigidez e vetores de cargas nodais equivalentes"
    Paint.dc:SetTextForeground(wx.wxBLACK)
    Paint.dc:DrawRotatedText("em coordenadas locais" , Paint.bar_canvas_end_x + 230, 1*146 + 8, 90)
    Paint.dc:DrawRotatedText("em coordenadas globais", Paint.bar_canvas_end_x + 230, 2*146 + 8, 90)
  end
  
  -- ALWAYS delete() any wxDCs created when done, to make sure it will reprint instead of printing over
  Paint.dc:delete()
end


--############################################################################################################################################################
-- WXLUA MAIN
--############################################################################################################################################################


-- Decimal to scientific notation converter (as a string) -- TODO: fix for numbers between -1 and 1
function Sci_notation(number, decimal_digits)
  local decimal = decimal_digits or 2
  local size = tonumber(string.len(tostring(number > 0 and number or -number)))
  
  local exponent = size - 1
  local significand = math.round(number/(10^exponent), 2)
  
  return fmt("%.2fE%d", significand, exponent)
end

-- Get correct angle in radians based on sin and cos
function Get_angle(sin, cos)
  local S, C = sin, cos
  local angle_rad, angle_deg
  
  if     S >= 0 and C >= 0 then angle_rad = math.asin(S)
  elseif S >= 0 and C < 0  then angle_rad = math.acos(C)
  elseif S < 0  and C <= 0 then angle_rad = math.pi - math.asin(S)
  elseif S <= 0 and C > 0  then angle_rad = 2*math.pi - math.acos(C)
  end
  
  angle_deg = 360*angle_rad/(2*math.pi)
  
  return angle_rad, angle_deg
end

-- Check if line is a point
function Is_point(xi, yi, xf, yf, round)
  if round then
    if math.round(xi) == math.round(xf) and math.round(yi) == math.round(yf) then return true else return false end
  else
    if xi == xf and yi == yf then return true else return false end
  end
end

-- Choice button to select bar
function DESLOCA.bar_choice_update()
  bar_choice_table = {} -- need to reset every time
  for i = 1, #DESLOCA.bar do
    bar_choice_table[i] = "Barra " .. i
  end
  
  if bar_choice ~= nil then
    bar_choice:Destroy()
  end
  bar_choice = wx.wxChoice(panel, IDs.choice_bar, wx.wxPoint(Paint.bar_canvas_x, 1), wx.wxSize(70, 20), bar_choice_table)
  
  -- Forcing a choice, to avoid it being blank at the start or when bar is created/deleted
  bar_choice:SetSelection(Paint.bar_sel-1)
end
DESLOCA.bar_choice_update()

-- Create/Delete bar buttons
local x_temp, y_temp = Paint.bar_canvas_x + bar_choice:GetSize():GetWidth() + 1, 0
create_bar_button = wx.wxButton(panel, IDs.create_bar_button, "Criar barra", wx.wxPoint(x_temp, y_temp))
x_temp = x_temp + create_bar_button:GetSize():GetWidth()
delete_bar_button = wx.wxButton(panel, IDs.delete_bar_button, "Excluir barra", wx.wxPoint(x_temp, y_temp))


-- Bar input/edit
local x_input, y_input = Paint.bar_canvas_x, Paint.bar_canvas_end_y + 8
x_temp, y_temp = x_input, y_input

local L, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
DESLOCA.bar[Paint.bar_sel].L = L -- since this was stored later only with events
local angle_rad, angle_deg = Get_angle(S, C)
local xi, yi, xf, yf = DESLOCA.bar[Paint.bar_sel].xi, DESLOCA.bar[Paint.bar_sel].yi, DESLOCA.bar[Paint.bar_sel].xf, DESLOCA.bar[Paint.bar_sel].yf
local delta_xi, delta_yi = xi, yi -- pra colocar bar em coordenadas locais
xi, yi, xf, yf = xi - delta_xi, yi - delta_yi, xf - delta_xi, yf - delta_yi -- pra colocar bar em coordenadas locais

static_text_L = wx.wxStaticText(panel, IDs.static_text_L, "L", wx.wxPoint(x_temp, y_temp), wx.wxSize(10, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 15
text_ctrl_L = wx.wxTextCtrl(panel, IDs.text_ctrl_L, tostring(L), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(50, 18)) ; x_temp = x_temp + 55
static_text_units_L = wx.wxStaticText(panel, IDs.static_text_units_L, "cm", wx.wxPoint(x_temp, y_temp), wx.wxSize(25, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 25

check_box_L_fixed = wx.wxCheckBox(panel, IDs.check_box_L_fixed, "fixo", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15)) ; x_temp, y_temp = x_input, y_temp + 20
--check_box_L_fixed:SetValue(false)

static_text_xiyi = wx.wxStaticText(panel, IDs.static_text_xiyi, "(xi, yi)", wx.wxPoint(x_temp, y_temp), wx.wxSize(35, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 40
text_ctrl_xi = wx.wxTextCtrl(panel, IDs.text_ctrl_xi, tostring(DESLOCA.bar[Paint.bar_sel].xi), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 40
text_ctrl_yi = wx.wxTextCtrl(panel, IDs.text_ctrl_yi, tostring(DESLOCA.bar[Paint.bar_sel].yi), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
static_text_units_xiyi = wx.wxStaticText(panel, IDs.static_text_units_xiyi, "cm", wx.wxPoint(x_temp, y_temp), wx.wxSize(25, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

static_text_xfxf = wx.wxStaticText(panel, IDs.static_text_xfxf, "(xf, yf)", wx.wxPoint(x_temp, y_temp), wx.wxSize(35, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 40
text_ctrl_xf = wx.wxTextCtrl(panel, IDs.text_ctrl_xf, tostring(DESLOCA.bar[Paint.bar_sel].xf), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 40
text_ctrl_yf = wx.wxTextCtrl(panel, IDs.text_ctrl_yf, tostring(DESLOCA.bar[Paint.bar_sel].yf), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
static_text_units_xfxf = wx.wxStaticText(panel, IDs.static_text_units_xfxf, "cm", wx.wxPoint(x_temp, y_temp), wx.wxSize(25, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

static_text_angle = wx.wxStaticText(panel, IDs.static_text_angle, "<", wx.wxPoint(x_temp, y_temp), wx.wxSize(10, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 15
spin_ctrl_angle = wx.wxSpinCtrl(panel, IDs.spin_ctrl_angle, tostring(angle_deg), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(50, 18), wx.wxSP_ARROW_KEYS, -1, 360) ; x_temp = x_temp + 55
static_text_units_angle = wx.wxStaticText(panel, IDs.static_text_units_angle, "°", wx.wxPoint(x_temp, y_temp), wx.wxSize(20, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input + 152, y_input

static_text_E = wx.wxStaticText(panel, IDs.static_text_E, "E", wx.wxPoint(x_temp, y_temp), wx.wxSize(10, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 15
text_ctrl_E = wx.wxTextCtrl(panel, IDs.text_ctrl_E, tostring(DESLOCA.bar[Paint.bar_sel].E), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(50, 18)) ; x_temp = x_temp + 55
static_text_units_E = wx.wxStaticText(panel, IDs.static_text_units_E, "kN/cm²", wx.wxPoint(x_temp, y_temp), wx.wxSize(48, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input + 152, y_temp + 20

static_text_A = wx.wxStaticText(panel, IDs.static_text_A, "A", wx.wxPoint(x_temp, y_temp), wx.wxSize(10, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 15
text_ctrl_A = wx.wxTextCtrl(panel, IDs.text_ctrl_A, tostring(DESLOCA.bar[Paint.bar_sel].A), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(50, 18)) ; x_temp = x_temp + 55
static_text_units_A = wx.wxStaticText(panel, IDs.static_text_units_A, "cm²", wx.wxPoint(x_temp, y_temp), wx.wxSize(30, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input + 152, y_temp + 20

static_text_I = wx.wxStaticText(panel, IDs.static_text_I, "I", wx.wxPoint(x_temp, y_temp), wx.wxSize(10, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 15
text_ctrl_I = wx.wxTextCtrl(panel, IDs.text_ctrl_I, tostring(DESLOCA.bar[Paint.bar_sel].I), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(50, 18)) ; x_temp = x_temp + 55
static_text_units_I = wx.wxStaticText(panel, IDs.static_text_units_I, "cm^4", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxALIGN_LEFT)


-- Distributed load input/edit
x_input, y_input = Paint.bar_canvas_end_x + 16, 22
x_temp, y_temp = x_input, y_input

wx.wxStaticText(panel, wx.wxID_ANY, "Carga distribuída", wx.wxPoint(x_temp, y_temp), wx.wxSize(90, 15), wx.wxALIGN_LEFT) ; y_temp = y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "q", wx.wxPoint(x_temp, y_temp), wx.wxSize(10, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 15
dist_load_text_ctrl_q = wx.wxTextCtrl(panel, IDs.dist_load_text_ctrl_q, tostring(DESLOCA.bar[Paint.bar_sel].dist.q), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
wx.wxStaticText(panel, wx.wxID_ANY, "kN/cm", wx.wxPoint(x_temp, y_temp), wx.wxSize(48, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

dist_load_option_local = wx.wxRadioButton(panel, IDs.dist_load_option_local, "local", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
dist_load_option_global_x = wx.wxRadioButton(panel, IDs.dist_load_option_global_x, "global em x", wx.wxPoint(x_temp, y_temp), wx.wxSize(75, 15)) ; x_temp = x_temp + 80
dist_load_option_global_y = wx.wxRadioButton(panel, IDs.dist_load_option_global_y, "global em y", wx.wxPoint(x_temp, y_temp), wx.wxSize(75, 15)) ; x_temp, y_temp = x_input, y_temp + 30
if DESLOCA.bar[Paint.bar_sel].dist.is_global == true then dist_load_option_global_x:SetValue(true) end


-- Nodal load input/edit
wx.wxStaticText(panel, wx.wxID_ANY, "Carga nodal", wx.wxPoint(x_temp, y_temp), wx.wxSize(70, 15), wx.wxALIGN_LEFT) ; y_temp = y_temp + 20

 -- for i node
nodal_load_static_text_i = wx.wxStaticText(panel, IDs.nodal_load_static_text_i, "Nó i (1)", wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, wx.wxALIGN_LEFT) ; y_temp = y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "Fx", wx.wxPoint(x_temp, y_temp), wx.wxSize(15, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 20
nodal_load_text_ctrl_Fx_i = wx.wxTextCtrl(panel, IDs.nodal_load_text_ctrl_Fx_i, tostring(DESLOCA.bar[Paint.bar_sel].i.Fx), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
wx.wxStaticText(panel, wx.wxID_ANY, "kN", wx.wxPoint(x_temp, y_temp), wx.wxSize(20, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "Fy", wx.wxPoint(x_temp, y_temp), wx.wxSize(15, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 20
nodal_load_text_ctrl_Fy_i = wx.wxTextCtrl(panel, IDs.nodal_load_text_ctrl_Fy_i, tostring(DESLOCA.bar[Paint.bar_sel].i.Fy), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
wx.wxStaticText(panel, wx.wxID_ANY, "kN", wx.wxPoint(x_temp, y_temp), wx.wxSize(20, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "Mz", wx.wxPoint(x_temp, y_temp), wx.wxSize(15, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 20
nodal_load_text_ctrl_Mz_i = wx.wxTextCtrl(panel, IDs.nodal_load_text_ctrl_Mz_i, tostring(DESLOCA.bar[Paint.bar_sel].i.Mz), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
wx.wxStaticText(panel, wx.wxID_ANY, "kNcm", wx.wxPoint(x_temp, y_temp), wx.wxSize(30, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

nodal_load_option_local_i = wx.wxRadioButton(panel, IDs.nodal_load_option_local_i, "local", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
nodal_load_option_global_i = wx.wxRadioButton(panel, IDs.nodal_load_option_global_i, "global", wx.wxPoint(x_temp, y_temp), wx.wxSize(50, 15)) ; x_temp, y_temp = x_input, y_temp + 20
if DESLOCA.bar[Paint.bar_sel].i.is_global == true then nodal_load_option_global_i:SetValue(true) end

 -- for f node
nodal_load_static_text_f = wx.wxStaticText(panel, IDs.nodal_load_static_text_f, "Nó f (2)", wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, wx.wxALIGN_LEFT) ; y_temp = y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "Fx", wx.wxPoint(x_temp, y_temp), wx.wxSize(15, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 20
nodal_load_text_ctrl_Fx_f = wx.wxTextCtrl(panel, IDs.nodal_load_text_ctrl_Fx_f, tostring(DESLOCA.bar[Paint.bar_sel].f.Fx), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
wx.wxStaticText(panel, wx.wxID_ANY, "kN", wx.wxPoint(x_temp, y_temp), wx.wxSize(20, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "Fy", wx.wxPoint(x_temp, y_temp), wx.wxSize(15, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 20
nodal_load_text_ctrl_Fy_f = wx.wxTextCtrl(panel, IDs.nodal_load_text_ctrl_Fy_f, tostring(DESLOCA.bar[Paint.bar_sel].f.Fy), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
wx.wxStaticText(panel, wx.wxID_ANY, "kN", wx.wxPoint(x_temp, y_temp), wx.wxSize(20, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "Mz", wx.wxPoint(x_temp, y_temp), wx.wxSize(15, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 20
nodal_load_text_ctrl_Mz_f = wx.wxTextCtrl(panel, IDs.nodal_load_text_ctrl_Mz_f, tostring(DESLOCA.bar[Paint.bar_sel].f.Mz), wx.wxPoint(x_temp, y_temp - 2), wx.wxSize(40, 18)) ; x_temp = x_temp + 45
wx.wxStaticText(panel, wx.wxID_ANY, "kNcm", wx.wxPoint(x_temp, y_temp), wx.wxSize(30, 15), wx.wxALIGN_LEFT) ; x_temp, y_temp = x_input, y_temp + 20

nodal_load_option_local_f = wx.wxRadioButton(panel, IDs.nodal_load_option_local_f, "local", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
nodal_load_option_global_f = wx.wxRadioButton(panel, IDs.nodal_load_option_global_f, "global", wx.wxPoint(x_temp, y_temp), wx.wxSize(50, 15)) ; x_temp, y_temp = x_input, y_temp + 20
if DESLOCA.bar[Paint.bar_sel].f.is_global == true then nodal_load_option_global_f:SetValue(true) end


-- Restraints input/edit
x_input, y_input = Paint.bar_canvas_end_x + 124, 22 + 70
x_temp, y_temp = x_input, y_input

wx.wxStaticText(panel, wx.wxID_ANY, "Vinculações", wx.wxPoint(x_temp, y_temp), wx.wxSize(60, 15), wx.wxALIGN_LEFT) ; y_temp = y_temp + 20

 -- for i node
restrain_static_text_i = wx.wxStaticText(panel, IDs.restrain_static_text_i, "Nó i (1)", wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, wx.wxALIGN_LEFT) ; y_temp = y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "X", wx.wxPoint(x_temp, y_temp), wx.wxSize(8, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 13
restrain_option_free_x_i = wx.wxRadioButton(panel, IDs.restrain_option_free_x_i, "livre", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
restrain_option_fixed_x_i = wx.wxRadioButton(panel, IDs.restrain_option_fixed_x_i, "fixo", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15)) ; x_temp, y_temp = x_input, y_temp + 20
if DESLOCA.bar[Paint.bar_sel].i.free_x == false then restrain_option_fixed_x_i:SetValue(true) end

wx.wxStaticText(panel, wx.wxID_ANY, "Y", wx.wxPoint(x_temp, y_temp), wx.wxSize(8, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 13
restrain_option_free_y_i = wx.wxRadioButton(panel, IDs.restrain_option_free_y_i, "livre", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
restrain_option_fixed_y_i = wx.wxRadioButton(panel, IDs.restrain_option_fixed_y_i, "fixo", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15)) ; x_temp, y_temp = x_input, y_temp + 20
if DESLOCA.bar[Paint.bar_sel].i.free_y == false then restrain_option_fixed_y_i:SetValue(true) end

wx.wxStaticText(panel, wx.wxID_ANY, "Z", wx.wxPoint(x_temp, y_temp), wx.wxSize(8, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 13
restrain_option_free_z_i = wx.wxRadioButton(panel, IDs.restrain_option_free_z_i, "livre", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
restrain_option_fixed_z_i = wx.wxRadioButton(panel, IDs.restrain_option_fixed_z_i, "fixo", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15)) ; x_temp, y_temp = x_input, y_temp + 40
if DESLOCA.bar[Paint.bar_sel].i.free_z == false then restrain_option_fixed_z_i:SetValue(true) end

 -- for i node
restrain_static_text_f = wx.wxStaticText(panel, IDs.restrain_static_text_f, "Nó f (2)", wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, wx.wxALIGN_LEFT) ; y_temp = y_temp + 20

wx.wxStaticText(panel, wx.wxID_ANY, "X", wx.wxPoint(x_temp, y_temp), wx.wxSize(8, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 13
restrain_option_free_x_f = wx.wxRadioButton(panel, IDs.restrain_option_free_x_f, "livre", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
restrain_option_fixed_x_f = wx.wxRadioButton(panel, IDs.restrain_option_fixed_x_f, "fixo", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15)) ; x_temp, y_temp = x_input, y_temp + 20
if DESLOCA.bar[Paint.bar_sel].f.free_x == false then restrain_option_fixed_x_f:SetValue(true) end

wx.wxStaticText(panel, wx.wxID_ANY, "Y", wx.wxPoint(x_temp, y_temp), wx.wxSize(8, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 13
restrain_option_free_y_f = wx.wxRadioButton(panel, IDs.restrain_option_free_y_f, "livre", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
restrain_option_fixed_y_f = wx.wxRadioButton(panel, IDs.restrain_option_fixed_y_f, "fixo", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15)) ; x_temp, y_temp = x_input, y_temp + 20
if DESLOCA.bar[Paint.bar_sel].f.free_y == false then restrain_option_fixed_y_f:SetValue(true) end

wx.wxStaticText(panel, wx.wxID_ANY, "Z", wx.wxPoint(x_temp, y_temp), wx.wxSize(8, 15), wx.wxALIGN_LEFT) ; x_temp = x_temp + 13
restrain_option_free_z_f = wx.wxRadioButton(panel, IDs.restrain_option_free_z_f, "livre", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15), wx.wxRB_GROUP) ; x_temp = x_temp + 45
restrain_option_fixed_z_f = wx.wxRadioButton(panel, IDs.restrain_option_fixed_z_f, "fixo", wx.wxPoint(x_temp, y_temp), wx.wxSize(40, 15)) ; x_temp, y_temp = x_input, y_temp + 20
if DESLOCA.bar[Paint.bar_sel].f.free_z == false then restrain_option_fixed_z_f:SetValue(true) end

-- Choice button for the bar display mode
x_temp, y_temp = Paint.bar_canvas_end_x + 230, 1
bar_mode_choice_table = {
  "Matrizes de rigidez e vetores de cargas nodais equivalentes",
  "Forças nas extremidades em coordenadas locais",
}
bar_mode_choice = wx.wxChoice(panel, IDs.bar_mode_choice, wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, bar_mode_choice_table)
bar_mode_choice:SetSelection(0) -- forcing a choice, to avoid it being blank at the start


-- Function that displays the bar stiffness matrixes, its load vectors and its rotation matrixes, if bar_mode_choice is "Matrizes de rigidez e vetores de cargas nodais equivalentes"
function DESLOCA.stiff_matrix_mode()
  
  -- Destroy controls from the other mode -- TODO: see if :Show(false) is better
  if static_text_f_ep ~= nil then static_text_f_ep:Destroy() end
  if grid_f_ep ~= nil then grid_f_ep:Destroy() end
  if static_text_mat_local ~= nil then static_text_mat_local:Destroy() end
  if grid_mat_local ~= nil then grid_mat_local:Destroy() end
  if static_text_R ~= nil then static_text_R:Destroy() end
  if grid_rotation_mat ~= nil then grid_rotation_mat:Destroy() end
  if static_text_u ~= nil then static_text_u:Destroy() end
  if grid_u ~= nil then grid_u:Destroy() end
  if static_text_f_local ~= nil then static_text_f_local:Destroy() end
  if grid_f_local ~= nil then grid_f_local:Destroy() end

  -- Bar stiffness matrix in local coordinates
  x_temp, y_temp = Paint.bar_canvas_end_x + 270, 30
  static_text_mat_local = wx.wxStaticText(panel, wx.wxID_ANY, "Matriz de Rigidez", wx.wxPoint(x_temp + 110, y_temp))
  static_text_mat_local_i_1 = wx.wxStaticText(panel, wx.wxID_ANY, "i", wx.wxPoint(x_temp - 15, y_temp + 38))
  static_text_mat_local_f_1 = wx.wxStaticText(panel, wx.wxID_ANY, "f", wx.wxPoint(x_temp - 15, y_temp + 89))
  static_text_mat_local_i_2 = wx.wxStaticText(panel, wx.wxID_ANY, "i", wx.wxPoint(x_temp + 73, y_temp + 5))
  static_text_mat_local_f_2 = wx.wxStaticText(panel, wx.wxID_ANY, "f", wx.wxPoint(x_temp + 222, y_temp + 5)) ; y_temp = y_temp + 20
  grid_mat_local = wx.wxGrid(panel, IDs.grid_mat_local, wx.wxPoint(x_temp, y_temp), wx.wxSize(315, 118), wx.wxALWAYS_SHOW_SB)
  grid_mat_local:CreateGrid(6, 6)
  grid_mat_local:EnableEditing(false)
  grid_mat_local:EnableDragGridSize(false)
  grid_mat_local:SetRowLabelSize(0)
  grid_mat_local:SetColLabelSize(0)
  grid_mat_local:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_mat_local:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_mat_local:SetDefaultColSize(50)
  grid_mat_local:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE)
  function DESLOCA.bar_local_grid_fill()
    for i = 0, 5 do
      for j = 0, 5 do
        local cell_value = fmt("%.1f", math.round(DESLOCA.bar[Paint.bar_sel].stiff_matrix_local[i+1][j+1], 1))
        grid_mat_local:SetCellValue(i, j, cell_value)
        if (i <= 2 and j >= 3) or (i >= 3 and j <= 2) then
          grid_mat_local:SetCellBackgroundColour(i, j, wx.wxColour(0xE0, 0xE0, 0xE0)) -- light gray
        end
      end
    end
  end

  -- Bar load vector in local coordinates
  x_temp, y_temp = x_temp + 292, y_temp - 20
  static_text_load_vector_local = wx.wxStaticText(panel, wx.wxID_ANY, "Vetor de cargas\n nodais equiv.", wx.wxPoint(x_temp + 12, y_temp - 10), wx.wxDefaultSize, wx.wxALIGN_CENTRE) ; y_temp = y_temp + 20
  grid_load_vector_local = wx.wxGrid(panel, IDs.grid_load_vector_local, wx.wxPoint(x_temp + 25, y_temp), wx.wxSize(60, 118), wx.wxALWAYS_SHOW_SB)
  grid_load_vector_local:CreateGrid(6, 1)
  grid_load_vector_local:EnableEditing(false)
  grid_load_vector_local:EnableDragGridSize(false)
  grid_load_vector_local:SetRowLabelSize(0)
  grid_load_vector_local:SetColLabelSize(0)
  grid_load_vector_local:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_load_vector_local:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_load_vector_local:SetDefaultColSize(50)
  --grid_load_vector_local:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE) -- TODO: decide
  --wx.wxStaticText(panel, wx.wxID_ANY, "", wx.wxPoint(x_temp + 70, y_temp), wx.wxSize(15, 100)) -- to hide disabled scroll arrows
  function DESLOCA.load_vector_local_grid_fill()
    for i = 0, 5 do
      local cell_value = fmt("%.1f", math.round(DESLOCA.bar[Paint.bar_sel].load_vector_local[i+1], 1))
      grid_load_vector_local:SetCellValue(i, 0, cell_value)
    end
  end

  -- Bar rotation matrixes
  x_temp, y_temp = x_temp + 90, y_temp - 20
  grid_rotation_mat_option_normal = wx.wxRadioButton(panel, IDs.grid_rotation_mat_option_normal, "Matriz de rotação", wx.wxPoint(x_temp + 12, y_temp - 20), wx.wxDefaultSize, wx.wxRB_GROUP)
  grid_rotation_mat_option_transposed = wx.wxRadioButton(panel, IDs.grid_rotation_mat_option_transposed, "Matriz de rotação transposta", wx.wxPoint(x_temp + 12, y_temp - 00), wx.wxDefaultSize) ; y_temp = y_temp + 20
  grid_rotation_mat = wx.wxGrid(panel, IDs.grid_rotation_mat, wx.wxPoint(x_temp, y_temp), wx.wxSize(196, 118), wx.wxALWAYS_SHOW_SB)
  grid_rotation_mat:CreateGrid(6, 6)
  grid_rotation_mat:EnableEditing(false)
  grid_rotation_mat:EnableDragGridSize(false)
  grid_rotation_mat:SetRowLabelSize(0)
  grid_rotation_mat:SetColLabelSize(0)
  grid_rotation_mat:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_rotation_mat:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_rotation_mat:SetDefaultColSize(30)
  grid_rotation_mat:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE)
  function DESLOCA.bar_rotation_grid_fill()
    for i = 0, 5 do
      for j = 0, 5 do
        local cell_value
        if grid_rotation_mat_option_normal:GetValue() == true then -- use grid for normal rotation matrix
          cell_value = fmt("%.2f", math.round(DESLOCA.bar[Paint.bar_sel].rotation_matrix[i+1][j+1], 2))
        else
          cell_value = fmt("%.2f", math.round(DESLOCA.bar[Paint.bar_sel].rotation_matrix_t[i+1][j+1], 2))
        end
        
        grid_rotation_mat:SetCellValue(i, j, cell_value)
      end
    end
  end

  -- Bar stiffness matrix in global coordinates
  x_temp, y_temp = Paint.bar_canvas_end_x + 270, y_temp + 18*7
  static_text_mat_global = wx.wxStaticText(panel, wx.wxID_ANY, "Matriz de Rigidez", wx.wxPoint(x_temp + 110, y_temp))
  node_info_static_text_i_1 = wx.wxStaticText(panel, IDs.node_info_static_text_i_1, "1", wx.wxPoint(x_temp - 15, y_temp + 38))
  node_info_static_text_f_1 = wx.wxStaticText(panel, IDs.node_info_static_text_f_1, "2", wx.wxPoint(x_temp - 15, y_temp + 89))
  node_info_static_text_i_2 = wx.wxStaticText(panel, IDs.node_info_static_text_i_2, "1", wx.wxPoint(x_temp + 70, y_temp + 5))
  node_info_static_text_f_2 = wx.wxStaticText(panel, IDs.node_info_static_text_f_2, "2", wx.wxPoint(x_temp + 222, y_temp + 5)) ; y_temp = y_temp + 20
  grid_mat_global = wx.wxGrid(panel, IDs.grid_mat_global, wx.wxPoint(x_temp, y_temp), wx.wxSize(315, 118), wx.wxALWAYS_SHOW_SB)
  grid_mat_global:CreateGrid(6, 6)
  grid_mat_global:EnableEditing(false)
  grid_mat_global:EnableDragGridSize(false)
  grid_mat_global:SetRowLabelSize(0)
  grid_mat_global:SetColLabelSize(0)
  grid_mat_global:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_mat_global:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_mat_global:SetDefaultColSize(50)
  grid_mat_global:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE)
  function DESLOCA.bar_global_grid_fill()
    for i = 0, 5 do
      for j = 0, 5 do
        local cell_value = fmt("%.1f", math.round(DESLOCA.bar[Paint.bar_sel].stiff_matrix_global[i+1][j+1], 1))
        grid_mat_global:SetCellValue(i, j, cell_value)
        if (i <= 2 and j >= 3) or (i >= 3 and j <= 2) then
          grid_mat_global:SetCellBackgroundColour(i, j, wx.wxColour(0xE0, 0xE0, 0xE0)) -- light gray
        end
      end
    end
  end

  -- Bar load vector in global coordinates
  x_temp, y_temp = x_temp + 292, y_temp -20
  static_text_load_vector_global = wx.wxStaticText(panel, wx.wxID_ANY, "Vetor de cargas\n nodais equiv.", wx.wxPoint(x_temp + 12, y_temp - 10), wx.wxDefaultSize, wx.wxALIGN_CENTRE) ; y_temp = y_temp + 20
  grid_load_vector_global = wx.wxGrid(panel, IDs.grid_load_vector_global, wx.wxPoint(x_temp + 25, y_temp), wx.wxSize(60, 118), wx.wxALWAYS_SHOW_SB)
  grid_load_vector_global:CreateGrid(6, 1)
  grid_load_vector_global:EnableEditing(false)
  grid_load_vector_global:EnableDragGridSize(false)
  grid_load_vector_global:SetRowLabelSize(0)
  grid_load_vector_global:SetColLabelSize(0)
  grid_load_vector_global:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_load_vector_global:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_load_vector_global:SetDefaultColSize(50)
  --grid_load_vector_global:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE) -- TODO: decide
  --wx.wxStaticText(panel, wx.wxID_ANY, "", wx.wxPoint(x_temp + 70, y_temp), wx.wxSize(15, 100)) -- to hide disabled scroll arrows
  function DESLOCA.load_vector_global_grid_fill()
    for i = 0, 5 do
      local cell_value = fmt("%.1f", math.round(DESLOCA.bar[Paint.bar_sel].load_vector_global[i+1], 1))
      grid_load_vector_global:SetCellValue(i, 0, cell_value)
    end
  end
  
end

-- Function that displays the matrixes and vectors used to get the forces on the ends of the bar, in an equation mode
function DESLOCA.forces_on_ends_mode()
  
  -- Destroy controls from the other mode -- TODO: see if :Show(false) is better
  grid_mat_local:Destroy()
  static_text_mat_local:Destroy()
  static_text_mat_local_i_1:Destroy()
  static_text_mat_local_f_1:Destroy()
  static_text_mat_local_i_2:Destroy()
  static_text_mat_local_f_2:Destroy()
  grid_load_vector_local:Destroy()
  static_text_load_vector_local:Destroy()
  grid_rotation_mat:Destroy()
  grid_rotation_mat_option_normal:Destroy()
  grid_rotation_mat_option_transposed:Destroy()
  grid_mat_global:Destroy()
  static_text_mat_global:Destroy()
  node_info_static_text_i_1:Destroy()
  node_info_static_text_f_1:Destroy()
  node_info_static_text_i_2:Destroy()
  node_info_static_text_f_2:Destroy()
  grid_load_vector_global:Destroy()
  static_text_load_vector_global:Destroy()
  
  -- Perfect settling forces vector in local coordinates
  x_temp, y_temp = Paint.bar_canvas_end_x + 230, 40
  static_text_f_ep = wx.wxStaticText(panel, wx.wxID_ANY, "f_ep\nlocal", wx.wxPoint(x_temp + 11, y_temp - 8), wx.wxDefaultSize, wx.wxALIGN_CENTRE) ; y_temp = y_temp + 20
  grid_f_ep = wx.wxGrid(panel, IDs.grid_f_ep, wx.wxPoint(x_temp, y_temp), wx.wxSize(60, 118), wx.wxALWAYS_SHOW_SB)
  grid_f_ep:CreateGrid(6, 1)
  grid_f_ep:EnableEditing(false)
  grid_f_ep:EnableDragGridSize(false)
  grid_f_ep:SetRowLabelSize(0)
  grid_f_ep:SetColLabelSize(0)
  grid_f_ep:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_f_ep:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_f_ep:SetDefaultColSize(50)
  function DESLOCA.f_ep_grid_fill()
    for i = 0, 5 do
      local cell_value = fmt("%.1f", math.round(-DESLOCA.bar[Paint.bar_sel].load_vector_local[i+1], 1))
      grid_f_ep:SetCellValue(i, 0, cell_value)
    end
  end
  
  -- Bar stiffness matrix in local coordinates
  x_temp = x_temp + 60 + 18
  static_text_mat_local = wx.wxStaticText(panel, wx.wxID_ANY, "K  local", wx.wxPoint(x_temp + 133, y_temp - 20))
  grid_mat_local = wx.wxGrid(panel, IDs.grid_mat_local, wx.wxPoint(x_temp, y_temp), wx.wxSize(315, 118), wx.wxALWAYS_SHOW_SB)
  grid_mat_local:CreateGrid(6, 6)
  grid_mat_local:EnableEditing(false)
  grid_mat_local:EnableDragGridSize(false)
  grid_mat_local:SetRowLabelSize(0)
  grid_mat_local:SetColLabelSize(0)
  grid_mat_local:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_mat_local:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_mat_local:SetDefaultColSize(50)
  grid_mat_local:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE)
  function DESLOCA.bar_local_grid_fill()
    for i = 0, 5 do
      for j = 0, 5 do
        local cell_value = fmt("%.1f", math.round(DESLOCA.bar[Paint.bar_sel].stiff_matrix_local[i+1][j+1], 1))
        grid_mat_local:SetCellValue(i, j, cell_value)
      end
    end
  end 

  -- Bar rotation matrix
  x_temp = x_temp + 315 + 12
  static_text_R = wx.wxStaticText(panel, wx.wxID_ANY, "R", wx.wxPoint(x_temp + 87, y_temp - 20))
  grid_rotation_mat = wx.wxGrid(panel, IDs.grid_rotation_mat, wx.wxPoint(x_temp, y_temp), wx.wxSize(196, 118), wx.wxALWAYS_SHOW_SB)
  grid_rotation_mat:CreateGrid(6, 6)
  grid_rotation_mat:EnableEditing(false)
  grid_rotation_mat:EnableDragGridSize(false)
  grid_rotation_mat:SetRowLabelSize(0)
  grid_rotation_mat:SetColLabelSize(0)
  grid_rotation_mat:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)
  grid_rotation_mat:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)
  grid_rotation_mat:SetDefaultColSize(30)
  grid_rotation_mat:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE)
  function DESLOCA.bar_rotation_grid_fill()
    for i = 0, 5 do
      for j = 0, 5 do
        local cell_value = fmt("%.2f", math.round(DESLOCA.bar[Paint.bar_sel].rotation_matrix[i+1][j+1], 2))
        grid_rotation_mat:SetCellValue(i, j, cell_value)
      end
    end
  end
  
  -- Bar displacements vector in global coordinates
  x_temp = x_temp + 196 + 12
  static_text_u = wx.wxStaticText(panel, wx.wxID_ANY, "u  global", wx.wxPoint(x_temp + 2, y_temp - 20), wx.wxDefaultSize, wx.wxALIGN_CENTRE)
  grid_u = wx.wxGrid(panel, IDs.grid_u, wx.wxPoint(x_temp, y_temp), wx.wxSize(60, 118), wx.wxALWAYS_SHOW_SB)
  grid_u:CreateGrid(6, 1)
  grid_u:EnableEditing(false)
  grid_u:EnableDragGridSize(false)
  grid_u:SetRowLabelSize(0)
  grid_u:SetColLabelSize(0)
  grid_u:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)
  grid_u:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)
  grid_u:SetDefaultColSize(55)
  function DESLOCA.u_grid_fill()
    for i = 0, 5 do
      local cell_value = fmt("%.6f", math.round(DESLOCA.bar[Paint.bar_sel].displacements_vector[i+1], 6))
      grid_u:SetCellValue(i, 0, cell_value)
    end
  end
  
  -- Bar forces on ends in local coordinates
  x_temp = x_temp + 60 + 18
  static_text_f_local = wx.wxStaticText(panel, wx.wxID_ANY, "f  local", wx.wxPoint(x_temp + 15, y_temp - 20), wx.wxDefaultSize, wx.wxALIGN_CENTRE)
  grid_f_local = wx.wxGrid(panel, IDs.grid_f_local, wx.wxPoint(x_temp, y_temp), wx.wxSize(78, 118), wx.wxALWAYS_SHOW_SB)
  grid_f_local:CreateGrid(6, 1)
  grid_f_local:EnableEditing(false)
  grid_f_local:EnableDragGridSize(false)
  grid_f_local:SetRowLabelSize(0)
  grid_f_local:SetColLabelSize(0)
  grid_f_local:SetScrollbar(wx.wxVERTICAL, 0, 0, 0)	
  grid_f_local:SetScrollbar(wx.wxHORIZONTAL, 0, 0, -1)	
  grid_f_local:SetDefaultColSize(50)
  function DESLOCA.f_local_grid_fill()
    for i = 0, 5 do
      local unit
      if i%3 + 1 == 1 then
        unit = "kN"
      elseif i%3 + 1 == 2 then
        unit = "kN"
      elseif i%3 + 1 == 3 then
        unit = "kNcm"
      end
      local cell_value = fmt("%.2f %s", math.round(DESLOCA.bar[Paint.bar_sel].forces_on_ends[i+1], 2), unit)
      grid_f_local:SetCellValue(i, 0, cell_value)
    end
    grid_f_local:AutoSizeColumns(true)
  end
  
end


-- Choice button for the structure display mode
structure_mode_choice_table = {
  "Estrutura com cargas",
  "Reações",
  "Configuração deformada",
  "Diagramas de esforços internos",
}
structure_mode_choice = wx.wxChoice(panel, IDs.structure_mode_choice, wx.wxPoint(Paint.canvas_x, Paint.canvas_y - 25), wx.wxDefaultSize, structure_mode_choice_table)
structure_mode_choice:SetSelection(0) -- forcing a choice, to avoid it being blank at the start
DESLOCA.structure_mode = 0 -- global value used in events and paint

-- Deformed scale slider only in deformed config mode
x_temp, y_temp = Paint.canvas_x + structure_mode_choice:GetSize():GetWidth() + 16, Paint.canvas_y - 32
deformed_scale_static_text = wx.wxStaticText(panel, IDs.deformed_scale_static_text, "Escala", wx.wxPoint(x_temp, y_temp)) ; x_temp, y_temp = x_temp + deformed_scale_static_text:GetSize():GetWidth() + 10, y_temp - 1
deformed_scale_slider = wx.wxSlider(panel, IDs.deformed_scale_slider, 2, 0, 4, wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, wx.wxSL_HORIZONTAL + wx.wxSL_AUTOTICKS) 
deformed_scale_static_text:Show(false) -- since the default structure mode is "Estrutura com cargas"
deformed_scale_slider:Show(false) -- since the default structure mode is "Estrutura com cargas"

-- Diagrams radio buttons only in diagrams mode
x_temp, y_temp = Paint.canvas_x + structure_mode_choice:GetSize():GetWidth() + 6, Paint.canvas_y - 20
diagram_option_N = wx.wxRadioButton(panel, IDs.diagram_option_N, "N (kN)", wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, wx.wxRB_GROUP) ; x_temp = x_temp + diagram_option_N:GetSize():GetWidth() + 5
diagram_option_Q = wx.wxRadioButton(panel, IDs.diagram_option_Q, "Q (kN)", wx.wxPoint(x_temp, y_temp)) ; x_temp = x_temp + diagram_option_Q:GetSize():GetWidth() + 5
diagram_option_M = wx.wxRadioButton(panel, IDs.diagram_option_M, "M (kNcm)", wx.wxPoint(x_temp, y_temp))
diagram_option_N:Show(false) -- since the default structure mode is "Estrutura com cargas"
diagram_option_Q:Show(false) -- since the default structure mode is "Estrutura com cargas"
diagram_option_M:Show(false) -- since the default structure mode is "Estrutura com cargas"

-- Function that handles the structure mode controls, based on the selected option
function DESLOCA.structure_mode_controls()
  -- Update the structure mode choice selected
  DESLOCA.structure_mode = structure_mode_choice:GetCurrentSelection()
  
  if DESLOCA.structure_mode == 0 then -- "Estrutura com cargas"
  
    deformed_scale_static_text:Show(false)
    deformed_scale_slider:Show(false)
    
    diagram_option_N:Show(false)
    diagram_option_Q:Show(false)
    diagram_option_M:Show(false)
    
  elseif DESLOCA.structure_mode == 1 then -- "Reações"
  
    deformed_scale_static_text:Show(false)
    deformed_scale_slider:Show(false)
    
    diagram_option_N:Show(false)
    diagram_option_Q:Show(false)
    diagram_option_M:Show(false)
    
  elseif DESLOCA.structure_mode == 2 then -- "Configuração deformada"
    
    deformed_scale_static_text:Show(true)
    deformed_scale_slider:Show(true)
    
    diagram_option_N:Show(false)
    diagram_option_Q:Show(false)
    diagram_option_M:Show(false)
    
  elseif DESLOCA.structure_mode == 3 then -- "Diagramas de esforços internos"
    
    deformed_scale_static_text:Show(false)
    deformed_scale_slider:Show(false)
    
    diagram_option_N:Show(true)
    diagram_option_Q:Show(true)
    diagram_option_M:Show(true)
    
  end
end

-- Function that create a deformed scale slider only in deformed config mode -- TODO: REMOVE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
function DESLOCA.deformed_config_mode()
  
  --wxSlider(wxWindow* parent, wxWindowID id, int value , int minValue, int maxValue, const wxPoint& point = wxDefaultPosition, const wxSize& size = wxDefaultSize, long style = wxSL_HORIZONTAL, const wxValidator& validator = wxDefaultValidator, const wxString& name = "slider")
  --x_temp, y_temp = Paint.canvas_x + structure_mode_choice:GetSize():GetWidth(), Paint.canvas_y - 25
  --deformed_scale_slider = wx.wxSlider(panel, IDs.deformed_scale_slider, 2, 0, 4, wx.wxPoint(x_temp, y_temp), wx.wxDefaultSize, wx.wxSL_HORIZONTAL + wx.wxSL_AUTOTICKS) 
  -- its event is wx.wxEVT_SCROLL_THUMBTRACK

end

-- Function that create force type radio buttons only in force diagrams mode
function DESLOCA.force_diagrams_mode()
  -- TODO
end

-- Structure stiffness matrix in global coordinates
x_temp, y_temp = Paint.canvas_end_x + 4, Paint.canvas_y - 20
static_text_mat_estrut = wx.wxStaticText(panel, wx.wxID_ANY, "Matriz de Rigidez da estrutura", wx.wxPoint(x_temp + 228, y_temp))
mat_estrut_simp_checkbox = wx.wxCheckBox(panel, IDs.mat_estrut_simp_checkbox, "Aplicar condições de contorno", wx.wxPoint(Paint.canvas_end_x + 228 + 190, y_temp)) ; y_temp = y_temp + 20
grid_mat_estrut = wx.wxGrid(panel, IDs.grid_mat_estrut, wx.wxPoint(x_temp, y_temp), wx.wxSize(496 + 3*40, 239 + 3*17), wx.wxALWAYS_SHOW_SB) -- TODO: resolve these values, should be related to the matrix size
grid_mat_estrut:CreateGrid(36, 36)
grid_mat_estrut:EnableEditing(false)
grid_mat_estrut:EnableDragGridSize(false)
grid_mat_estrut:SetRowLabelSize(45)
grid_mat_estrut:SetColLabelSize(18)
grid_mat_estrut:SetDefaultColSize(45)
grid_mat_estrut:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE)
for j = 0, grid_mat_estrut:GetNumberCols() do -- change col label to numbers instead of default letters
  --[[
  if (j)%3 + 1 == 1 then
    grid_mat_estrut:SetColLabelValue(j, fmt("%d x", math.ceil((j+1)/3)))
    grid_mat_estrut:SetRowLabelValue(j, fmt("%d x", math.ceil((j+1)/3)))
  elseif (j)%3 + 1 == 2 then
    grid_mat_estrut:SetColLabelValue(j, fmt("%d y", math.ceil((j+1)/3)))
    grid_mat_estrut:SetRowLabelValue(j, fmt("%d y", math.ceil((j+1)/3)))
  elseif (j)%3 + 1 == 3 then
    grid_mat_estrut:SetColLabelValue(j, fmt("%d z", math.ceil((j+1)/3)))
    grid_mat_estrut:SetRowLabelValue(j, fmt("%d z", math.ceil((j+1)/3)))
  end
  ]]
  grid_mat_estrut:SetColLabelValue(j, fmt("%d (%d)", j+1, math.ceil((j+1)/3)))
  grid_mat_estrut:SetRowLabelValue(j, fmt("%d (%d)", j+1, math.ceil((j+1)/3)))
end
function DESLOCA.structure_grid_fill(format)
  for i = 0, 36 - 1 do
    for j = 0, 36 - 1 do
      -- Get cell value
      local value_str
      
      if mat_estrut_simp_checkbox:GetValue() == false then -- normal matrix -- TODO: save some lines here
        if i >= #DESLOCA.structure_stiff_matrix or j >= #DESLOCA.structure_stiff_matrix then -- Clear "out-of-bounds" cells
          value_str = ""
        else -- Get cell value
          if format == "scientific" then
            value_str = Sci_notation(DESLOCA.structure_stiff_matrix[i+1][j+1])
          else
            value_str = fmt("%.1f", math.round(DESLOCA.structure_stiff_matrix[i+1][j+1], 1))
          end
        end
      else -- simplified matrix
        if i >= #DESLOCA.structure_stiff_matrix_simp or j >= #DESLOCA.structure_stiff_matrix_simp then -- Clear "out-of-bounds" cells
          value_str = ""
        else -- Get cell value
          if format == "scientific" then
            value_str = Sci_notation(DESLOCA.structure_stiff_matrix_simp[i+1][j+1])
          else
            value_str = fmt("%.1f", math.round(DESLOCA.structure_stiff_matrix_simp[i+1][j+1], 1))
          end
        end
      end
      
      -- Handle indetermination
      if DESLOCA.indeterminate_system then
        value_str = ""
      end
      
      -- Write value
      grid_mat_estrut:SetCellValue(i, j, value_str)
      
      -- Set default cell formatting
      grid_mat_estrut:SetCellFont(i, j, wx.wxFont(7, wx.wxFONTFAMILY_DEFAULT, wx.wxFONTSTYLE_NORMAL, wx.wxNORMAL))
      grid_mat_estrut:SetCellBackgroundColour(i, j, wx.wxWHITE)
      grid_mat_estrut:SetCellTextColour(i, j, wx.wxBLACK)
      
      -- Paint checkered background
      if (math.ceil((i+1)/3)%2 == 1 and math.ceil((j+1)/3)%2 == 0) or (math.ceil((i+1)/3)%2 == 0 and math.ceil((j+1)/3)%2 == 1) then
        grid_mat_estrut:SetCellBackgroundColour(i, j, wx.wxColour(0xE0, 0xE0, 0xE0)) -- light gray
      end
      
      -- Paint selected bar cells
      local sel_node_i, sel_node_f = DESLOCA.bar[Paint.bar_sel].global_node.i, DESLOCA.bar[Paint.bar_sel].global_node.f
      local sel_colour = wx.wxColour("ORANGE")
      -- ii
      if i+1 >= 3*sel_node_i - 2 and i+1 <= 3*sel_node_i and j+1 >= 3*sel_node_i - 2 and j+1 <= 3*sel_node_i then
        grid_mat_estrut:SetCellBackgroundColour(i, j, sel_colour)
        grid_mat_estrut:SetCellTextColour(i, j, wx.wxWHITE)        
      end
      -- if
      if i+1 >= 3*sel_node_i - 2 and i+1 <= 3*sel_node_i and j+1 >= 3*sel_node_f - 2 and j+1 <= 3*sel_node_f then
        grid_mat_estrut:SetCellBackgroundColour(i, j, sel_colour)
        grid_mat_estrut:SetCellTextColour(i, j, wx.wxWHITE)        
      end
      -- fi
      if i+1 >= 3*sel_node_f - 2 and i+1 <= 3*sel_node_f and j+1 >= 3*sel_node_i - 2 and j+1 <= 3*sel_node_i then
        grid_mat_estrut:SetCellBackgroundColour(i, j, sel_colour)
        grid_mat_estrut:SetCellTextColour(i, j, wx.wxWHITE)        
      end
      -- ff
      if i+1 >= 3*sel_node_f - 2 and i+1 <= 3*sel_node_f and j+1 >= 3*sel_node_f - 2 and j+1 <= 3*sel_node_f then
        grid_mat_estrut:SetCellBackgroundColour(i, j, sel_colour)
        grid_mat_estrut:SetCellTextColour(i, j, wx.wxWHITE)        
      end
    end
  end
  
  -- Make selected bar cells as visible as possible
  local sel_node_i, sel_node_f = DESLOCA.bar[Paint.bar_sel].global_node.i, DESLOCA.bar[Paint.bar_sel].global_node.f
  local highest_node = sel_node_i > sel_node_f and sel_node_i or sel_node_f
  local lowest_node = sel_node_i < sel_node_f and sel_node_i or sel_node_f
  if (3*highest_node - 0) - (3*lowest_node - 2) < 15 then grid_mat_estrut:MakeCellVisible(3*highest_node - 0 - 1, 3*highest_node - 0 - 1) end
  grid_mat_estrut:MakeCellVisible(3*lowest_node - 2 - 1, 3*lowest_node - 2 - 1)
end

-- Structure load vector in global coordinates
x_temp, y_temp = x_temp + 496 + 3*40, y_temp + 18
static_text_load_vector_structure = wx.wxStaticText(panel, wx.wxID_ANY, "Vetor de cargas nodais\nequiv. da estrutura", wx.wxPoint(x_temp + 12, y_temp - 80), wx.wxDefaultSize, wx.wxALIGN_CENTRE)-- ; y_temp = y_temp + 20
grid_load_vector_structure_option1 = wx.wxRadioButton(panel, IDs.grid_load_vector_structure_option1, "Resultante", wx.wxPoint(x_temp + 4, y_temp - 50), wx.wxDefaultSize, wx.wxRB_GROUP)
grid_load_vector_structure_option2 = wx.wxRadioButton(panel, IDs.grid_load_vector_structure_option2, "Parcela cargas distr.", wx.wxPoint(x_temp + 4, y_temp - 35), wx.wxDefaultSize)
grid_load_vector_structure_option3 = wx.wxRadioButton(panel, IDs.grid_load_vector_structure_option3, "Parcela cargas nodais", wx.wxPoint(x_temp + 4, y_temp - 20), wx.wxDefaultSize)
grid_load_vector_structure = wx.wxGrid(panel, IDs.grid_load_vector_structure, wx.wxPoint(x_temp + 14, y_temp), wx.wxSize(106, 221 + 3*17), wx.wxALWAYS_SHOW_SB) -- TODO: resolve this "221", should be related to the matrix size
grid_load_vector_structure:CreateGrid(36, 1)
grid_load_vector_structure:EnableEditing(false)
grid_load_vector_structure:EnableDragGridSize(false)
grid_load_vector_structure:SetRowLabelSize(25)
grid_load_vector_structure:SetColLabelSize(0)
grid_load_vector_structure:SetDefaultColSize(66)
--grid_load_vector_structure:SetDefaultCellAlignment(wx.wxALIGN_CENTRE, wx.wxALIGN_CENTRE) -- TODO: decide
function DESLOCA.structure_load_vector_grid_fill()
  for i = 0, 36 - 1 do
    -- Get cell value
    local value_str
    local simplified = mat_estrut_simp_checkbox:GetValue()
    if i >= #DESLOCA.structure_stiff_matrix then -- clear "out-of-bounds" cells
      value_str = ""
    else -- get cell value
    
      if grid_load_vector_structure_option1:GetValue() == true then -- show vetor result
      
        value_str = fmt("%.1f", simplified and math.round(DESLOCA.structure_load_vector_simp[i+1], 1) or math.round(DESLOCA.structure_load_vector[i+1], 1))
    
      elseif grid_load_vector_structure_option2:GetValue() == true then -- show distributed load portion
      
        value_str = fmt("%.1f", math.round(DESLOCA.structure_load_vector_dist_portion[i+1], 1))
      
      elseif grid_load_vector_structure_option3:GetValue() == true then -- show nodal load portion
      
        value_str = fmt("%.1f", math.round(DESLOCA.structure_load_vector_nodal_portion[i+1], 1))
      
      end
    end
    
    -- Handle indetermination
    if DESLOCA.indeterminate_system then
      value_str = ""
    end
    
    -- Write value
    grid_load_vector_structure:SetCellValue(i, 0, value_str)
    
    -- Set default cell formatting
    grid_load_vector_structure:SetCellBackgroundColour(i, 0, wx.wxWHITE)
    grid_load_vector_structure:SetCellTextColour(i, 0, wx.wxBLACK)
    
    -- Paint checkered background
    if math.ceil((i+1)/3)%2 == 1 then
      grid_load_vector_structure:SetCellBackgroundColour(i, 0, wx.wxColour(0xE0, 0xE0, 0xE0)) -- light gray
    end
    
    -- Paint selected bar cells
    local sel_node_i, sel_node_f = DESLOCA.bar[Paint.bar_sel].global_node.i, DESLOCA.bar[Paint.bar_sel].global_node.f
    local sel_colour = wx.wxColour("ORANGE")
    -- Fi
    if i+1 >= 3*sel_node_i - 2 and i+1 <= 3*sel_node_i then
      grid_load_vector_structure:SetCellBackgroundColour(i, 0, sel_colour)
      grid_load_vector_structure:SetCellTextColour(i, 0, wx.wxWHITE)        
    end
    -- Ff
    if i+1 >= 3*sel_node_f - 2 and i+1 <= 3*sel_node_f then
      grid_load_vector_structure:SetCellBackgroundColour(i, 0, sel_colour)
      grid_load_vector_structure:SetCellTextColour(i, 0, wx.wxWHITE)        
    end
  end
  
  -- Make selected bar cells as visible as possible
  local sel_node_i, sel_node_f = DESLOCA.bar[Paint.bar_sel].global_node.i, DESLOCA.bar[Paint.bar_sel].global_node.f
  local highest_node = sel_node_i > sel_node_f and sel_node_i or sel_node_f
  local lowest_node = sel_node_i < sel_node_f and sel_node_i or sel_node_f
  if (3*highest_node - 0) - (3*lowest_node - 2) < 15 then grid_load_vector_structure:MakeCellVisible(3*highest_node - 0 - 1, 0) end
  grid_load_vector_structure:MakeCellVisible(3*lowest_node - 2 - 1, 0)
end

-- Displacements and Reactions
x_temp, y_temp = x_temp + 120, y_temp
grid_results_option_disp = wx.wxRadioButton(panel, IDs.grid_results_option_disp, "Deslocamentos", wx.wxPoint(x_temp + 50, y_temp - 40), wx.wxDefaultSize, wx.wxRB_GROUP)
grid_results_option_react = wx.wxRadioButton(panel, IDs.grid_results_option_react, "Reações", wx.wxPoint(x_temp + 50, y_temp - 20), wx.wxDefaultSize)
grid_results = wx.wxGrid(panel, IDs.grid_results, wx.wxPoint(x_temp + 14, y_temp), wx.wxSize(160, 221 + 3*17), wx.wxALWAYS_SHOW_SB)
grid_results:CreateGrid(36, 1)
grid_results:EnableEditing(false)
grid_results:EnableDragGridSize(false)
grid_results:SetRowLabelSize(25)
grid_results:SetColLabelSize(0)
grid_results:SetDefaultColSize(150)
grid_results:SetRowLabelSize(30)
function DESLOCA.results_grid_fill()
  for i = 0, 36 - 1 do
    -- Get cell value
    local value_str
    if i >= #DESLOCA.structure_stiff_matrix then -- clear "out-of-bounds" cells
      value_str = ""
      grid_results:SetRowLabelValue(i, "")
    else -- get cell value
    
      if grid_results_option_disp:GetValue() == true then -- use grid for displacements
      
        grid_results:SetRowLabelValue(i, fmt("D%d", i+1))
        
        local displacement_direction, unit
        if i%3 + 1 == 1 then 
          displacement_direction = "x"
          unit = "cm"
        elseif i%3 + 1 == 2 then
          displacement_direction = "y"
          unit = "cm"
        elseif i%3 + 1 == 3 then
          displacement_direction = "z"
          unit = "rad"
        end
        value_str = fmt(" %s%d = %.6f %s", displacement_direction, math.ceil((i+1)/3), math.round(DESLOCA.displacements[i+1], 6), unit)
        
      elseif grid_results_option_react:GetValue() == true then -- use grid for reactions
      
        grid_results:SetRowLabelValue(i, fmt("R%d", i+1))
        
        local reaction_direction, unit
        if DESLOCA.node[math.ceil((i+1)/3)].free_x == false and i%3 + 1 == 1 then 
          reaction_direction = "H"
          unit = "kN"
        elseif DESLOCA.node[math.ceil((i+1)/3)].free_y == false and i%3 + 1 == 2 then
          reaction_direction = "V"
          unit = "kN"
        elseif DESLOCA.node[math.ceil((i+1)/3)].free_z == false and i%3 + 1 == 3 then
          reaction_direction = "M"
          unit = "kNcm"
        end
        if reaction_direction ~= nil then -- has restraint, so reaction makes sense
          value_str = fmt(" %s%d = %.2f %s", reaction_direction, math.ceil((i+1)/3), math.round(DESLOCA.reactions[i+1], 2), unit)
        else -- hasn't, so no purpose of showing the value here
          value_str = " -------"
        end
        
      end
    end
    
    -- Handle indetermination
    if DESLOCA.indeterminate_system then
      value_str = ""
    end
    
    -- Write value
    grid_results:SetCellValue(i, 0, value_str)
    
    -- Set default cell formatting
    grid_results:SetCellBackgroundColour(i, 0, wx.wxWHITE)
    grid_results:SetCellTextColour(i, 0, wx.wxBLACK)
    
    -- Paint checkered background
    if math.ceil((i+1)/3)%2 == 1 then
      grid_results:SetCellBackgroundColour(i, 0, wx.wxColour(0xE0, 0xE0, 0xE0)) -- light gray
    end
    
    -- Paint selected bar cells
    local sel_node_i, sel_node_f = DESLOCA.bar[Paint.bar_sel].global_node.i, DESLOCA.bar[Paint.bar_sel].global_node.f
    local sel_colour = wx.wxColour("ORANGE")
    -- Fi
    if i+1 >= 3*sel_node_i - 2 and i+1 <= 3*sel_node_i then
      grid_results:SetCellBackgroundColour(i, 0, sel_colour)
      grid_results:SetCellTextColour(i, 0, wx.wxWHITE)        
    end
    -- Ff
    if i+1 >= 3*sel_node_f - 2 and i+1 <= 3*sel_node_f then
      grid_results:SetCellBackgroundColour(i, 0, sel_colour)
      grid_results:SetCellTextColour(i, 0, wx.wxWHITE)        
    end
  end
  
  -- Adjust colum size to fit the largest result string
  grid_results:AutoSizeColumns(true)
  
  -- Make selected bar cells as visible as possible -- TODO:
  local sel_node_i, sel_node_f = DESLOCA.bar[Paint.bar_sel].global_node.i, DESLOCA.bar[Paint.bar_sel].global_node.f
  local highest_node = sel_node_i > sel_node_f and sel_node_i or sel_node_f
  local lowest_node = sel_node_i < sel_node_f and sel_node_i or sel_node_f
  if (3*highest_node - 0) - (3*lowest_node - 2) < 15 then grid_results:MakeCellVisible(3*highest_node - 0 - 1, 0) end
  grid_results:MakeCellVisible(3*lowest_node - 2 - 1, 0)
end


--############################################################################################################################################################
-- EVENTS
--############################################################################################################################################################

--- FUNCTIONS

-- Function to update the whole program, used by all the events
function DESLOCA.update_all()
  -- Evaluate nodes
  DESLOCA.eval_nodes()
  
  -- Load deafult structure stiffness matrixes
  DESLOCA.all_matrix_update()
  
  -- Get default structure load vectors
  DESLOCA.load_vector_update()
  
  -- Get the displacements
  DESLOCA.get_displacements()
  
  -- Get the reactions
  DESLOCA.get_reactions()
  
  -- Get the forces on the bar ends
  DESLOCA.get_forces_on_ends()
  
  -- Get N, Q and M for the forces diagrams
  DESLOCA.get_diagram_forces()
  
  -- Fill grids with stiffness matrixes
  DESLOCA.grid_fill_all()
  
  -- Update the node info
  DESLOCA.node_info_update()
  
  -- Refresh the program window
  frame:Refresh()
  
  -- Clear status bar
  DESLOCA.clear_status_bar(true, true)
end

-- Function to update the bar info when you change between them or create/delete one -- MUST be called after update_all() due node_info_update()
function DESLOCA.update_bar_info()
  Paint.bar_sel = bar_choice:GetCurrentSelection() + 1
  
  --DESLOCA.grid_fill_all()
  
  -- TODO: USE ANOTHER FUNCTION FOR THESE, PROBABLY THIS ONE (resolve lines 324~326)
  local L, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
  DESLOCA.bar[Paint.bar_sel].L = L -- to make sure
  local angle_rad, angle_deg = Get_angle(S, C)
  local xi, yi, xf, yf = DESLOCA.bar[Paint.bar_sel].xi, DESLOCA.bar[Paint.bar_sel].yi, DESLOCA.bar[Paint.bar_sel].xf, DESLOCA.bar[Paint.bar_sel].yf
  
  -- Properties
  text_ctrl_L:ChangeValue(tostring(L))
  text_ctrl_xi:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xi))
  text_ctrl_yi:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yi))
  text_ctrl_xf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xf))
  text_ctrl_yf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yf))
  spin_ctrl_angle:SetValue(angle_deg)
  text_ctrl_E:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].E))
  text_ctrl_A:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].A))
  text_ctrl_I:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].I))
  
  -- Distributed load
  dist_load_text_ctrl_q:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].dist.q))
  if DESLOCA.bar[Paint.bar_sel].dist.is_global == false then
    dist_load_option_local:SetValue(true)
  else
    if DESLOCA.bar[Paint.bar_sel].dist.global_dir == "x" then
      dist_load_option_global_x:SetValue(true)
    else
      dist_load_option_global_y:SetValue(true)
    end
  end
  
  -- Nodal loads
  nodal_load_text_ctrl_Fx_i:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].i.Fx))
  nodal_load_text_ctrl_Fy_i:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].i.Fy))
  nodal_load_text_ctrl_Mz_i:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].i.Mz))
  if DESLOCA.bar[Paint.bar_sel].i.is_global == false then nodal_load_option_local_i:SetValue(true) else nodal_load_option_global_i:SetValue(true) end
  nodal_load_text_ctrl_Fx_f:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].f.Fx))
  nodal_load_text_ctrl_Fy_f:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].f.Fy))
  nodal_load_text_ctrl_Mz_f:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].f.Mz))
  if DESLOCA.bar[Paint.bar_sel].f.is_global == false then nodal_load_option_local_f:SetValue(true) else nodal_load_option_global_f:SetValue(true) end
  
  -- Restraints
  if DESLOCA.bar[Paint.bar_sel].i.free_x == false then restrain_option_fixed_x_i:SetValue(true) else restrain_option_free_x_i:SetValue(true) end
  if DESLOCA.bar[Paint.bar_sel].i.free_y == false then restrain_option_fixed_y_i:SetValue(true) else restrain_option_free_y_i:SetValue(true) end
  if DESLOCA.bar[Paint.bar_sel].i.free_z == false then restrain_option_fixed_z_i:SetValue(true) else restrain_option_free_z_i:SetValue(true) end
  if DESLOCA.bar[Paint.bar_sel].f.free_x == false then restrain_option_fixed_x_f:SetValue(true) else restrain_option_free_x_f:SetValue(true) end
  if DESLOCA.bar[Paint.bar_sel].f.free_y == false then restrain_option_fixed_y_f:SetValue(true) else restrain_option_free_y_f:SetValue(true) end
  if DESLOCA.bar[Paint.bar_sel].f.free_z == false then restrain_option_fixed_z_f:SetValue(true) else restrain_option_free_z_f:SetValue(true) end
  
  -- Node info
  DESLOCA.node_info_update()
  
  -- Refresh the drawings and the program itself
  frame:Refresh()
end

-- Function to fill all the grids used in the program
function DESLOCA.grid_fill_all()
  if bar_mode_choice:GetCurrentSelection() == 0 then -- "Matrizes de rigidez e vetores de cargas nodais equivalentes"
  
    DESLOCA.bar_local_grid_fill()
    DESLOCA.load_vector_local_grid_fill()
    
    DESLOCA.bar_rotation_grid_fill()
    
    DESLOCA.bar_global_grid_fill()
    DESLOCA.load_vector_global_grid_fill()
    
  elseif bar_mode_choice:GetCurrentSelection() == 1 then -- "Forças nas extremidades em coordenadas locais"
  
    DESLOCA.f_ep_grid_fill()
    DESLOCA.bar_local_grid_fill()
    DESLOCA.bar_rotation_grid_fill()
    DESLOCA.u_grid_fill()
    
    DESLOCA.f_local_grid_fill()
    
  end
  
  DESLOCA.structure_grid_fill()
  DESLOCA.structure_load_vector_grid_fill()
  
  DESLOCA.results_grid_fill()
end

-- Function that updates the node info in static texts
function DESLOCA.node_info_update()
  if bar_mode_choice:GetCurrentSelection() ~= 0 then return end
  
  nodal_load_static_text_i:SetLabel(fmt("Nó i (%d)", DESLOCA.bar[Paint.bar_sel].global_node.i))
  nodal_load_static_text_f:SetLabel(fmt("Nó f (%d)", DESLOCA.bar[Paint.bar_sel].global_node.f))
  
  restrain_static_text_i:SetLabel(fmt("Nó i (%d)", DESLOCA.bar[Paint.bar_sel].global_node.i))
  restrain_static_text_f:SetLabel(fmt("Nó f (%d)", DESLOCA.bar[Paint.bar_sel].global_node.f))
  
  node_info_static_text_i_1:SetLabel(fmt("%d", DESLOCA.bar[Paint.bar_sel].global_node.i))
  node_info_static_text_f_1:SetLabel(fmt("%d", DESLOCA.bar[Paint.bar_sel].global_node.f))
  
  node_info_static_text_i_2:SetLabel(fmt("%d", DESLOCA.bar[Paint.bar_sel].global_node.i))
  node_info_static_text_f_2:SetLabel(fmt("%d", DESLOCA.bar[Paint.bar_sel].global_node.f)) 
end

-- Function that get user cell click and display the whole cell value in the status bar
function DESLOCA.grid_click_handle(event, self)
  local sel_row, sel_col = event:GetRow(), event:GetCol()
  self:ClearSelection() -- to avoid persistance of block selection 
  self:SetGridCursor(sel_row, sel_col) -- to avoid selection always going to cell 0,0
  status_bar:SetStatusText(self:GetCellValue(sel_row, sel_col), 1)
  
  -- Export value to clipboard
  local clipBoard = wx.wxClipboard.Get()
  if clipBoard and clipBoard:Open() then
    clipBoard:SetData(wx.wxTextDataObject(self:GetCellValue(sel_row, sel_col)))
    clipBoard:Flush()
    clipBoard:Close()
  end
end

function DESLOCA.clear_status_bar(field0, field1, field2)
  if field0 == true then status_bar:SetStatusText("", 0) end
  if field1 == true then status_bar:SetStatusText("", 1) end
  if field2 == true then status_bar:SetStatusText("", 2) end
end

--- EVENTS

-- Connect the Paint main function with the paint event
panel:Connect(wx.wxEVT_PAINT, Paint.paint)

-- Bar choice button
panel:Connect(IDs.choice_bar, wx.wxEVT_COMMAND_CHOICE_SELECTED, function()
  DESLOCA.update_bar_info()
  DESLOCA.grid_fill_all()
  DESLOCA.clear_status_bar(true, true)
end)

-- Bar create button
panel:Connect(IDs.create_bar_button, wx.wxEVT_COMMAND_BUTTON_CLICKED, function()
  local structure_size = #DESLOCA.bar
  local bar_limit = 10
  if structure_size < bar_limit then
    DESLOCA.bar[structure_size + 1] = {xi = DESLOCA.bar[Paint.bar_sel].xf, yi = DESLOCA.bar[Paint.bar_sel].yf, xf = DESLOCA.bar[Paint.bar_sel].xf + 200, yf = DESLOCA.bar[Paint.bar_sel].yf,
      E = DESLOCA.bar[Paint.bar_sel].E, A = DESLOCA.bar[Paint.bar_sel].A, I = DESLOCA.bar[Paint.bar_sel].I, global_node = {}, dist = {q = 0, is_global = false, global_dir = "x"},
      i = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = true, free_y = true, free_z = true}, f = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = true, free_y = true, free_z = true}}
    Paint.bar_sel = structure_size + 1
    local L = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
    DESLOCA.bar[Paint.bar_sel].L = L -- to make sure
    DESLOCA.bar_choice_update()
    DESLOCA.update_all()
    DESLOCA.update_bar_info()
    frame:SetStatusText("Barra adicionada!")
  else
    wx.wxMessageBox(fmt("O DESLOCA aceita até %d barras.", bar_limit), "Aviso!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  end
end)

-- Bar delete button
panel:Connect(IDs.delete_bar_button, wx.wxEVT_COMMAND_BUTTON_CLICKED, function()
  local structure_size = #DESLOCA.bar
  if structure_size > 1 then
    for i = Paint.bar_sel, structure_size do
      DESLOCA.bar[i] = DESLOCA.bar[i + 1] -- TODO: check if this is safe, since there's no table copying in Lua
    end
    if Paint.bar_sel == structure_size then -- deleted last bar in the list, conditional needed to avoid reading nil
      Paint.bar_sel = structure_size - 1
    end
    DESLOCA.bar_choice_update()
    DESLOCA.update_all()
    DESLOCA.update_bar_info()
    frame:SetStatusText("Barra excluída!")
  else
    wx.wxMessageBox("Não é possível deixar a estrutura sem nenhuma barra.", "Aviso!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  end
end)

-- Bar edit "L"
panel:Connect(IDs.text_ctrl_L, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local L_new = tonumber(text_ctrl_L:GetValue()) -- accept only numbers
  if text_ctrl_L:GetValue() == "" or text_ctrl_L:GetValue() == "-" or text_ctrl_L:GetValue() == "." or text_ctrl_L:GetValue() == "-." then
    -- wait number or error
  elseif L_new == nil then
    wx.wxMessageBox("O valor inserido para \"L\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  elseif L_new < 0 then
    wx.wxMessageBox("Barra com comprimento \"negativo\", insira novo valor.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
    text_ctrl_L:ChangeValue("")
  elseif L_new < 1 then
    wx.wxMessageBox("Barra muito pequena, insira novo comprimento.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
    text_ctrl_L:ChangeValue("")
  else
    local _, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
    DESLOCA.bar[Paint.bar_sel].L = L_new
    DESLOCA.bar[Paint.bar_sel].xf = DESLOCA.bar[Paint.bar_sel].xi + L_new*C
    DESLOCA.bar[Paint.bar_sel].yf = DESLOCA.bar[Paint.bar_sel].yi + L_new*S
    text_ctrl_xf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xf))
    text_ctrl_yf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yf))
    DESLOCA.update_all()
  end
end)

-- Bar edit "xi"
panel:Connect(IDs.text_ctrl_xi, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local xi_new = tonumber(text_ctrl_xi:GetValue()) -- accept only numbers
  if text_ctrl_xi:GetValue() == "" or text_ctrl_xi:GetValue() == "-" or text_ctrl_xi:GetValue() == "." or text_ctrl_xi:GetValue() == "-." then
    -- wait number or error
  elseif xi_new == nil then
    wx.wxMessageBox("O valor inserido para \"xi\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  elseif Is_point(xi_new, DESLOCA.bar[Paint.bar_sel].yi, DESLOCA.bar[Paint.bar_sel].xf, DESLOCA.bar[Paint.bar_sel].yf, true) then
    wx.wxMessageBox("Barra muito pequena ou pontual, insira novas coordenadas.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
    text_ctrl_xi:ChangeValue("")
  else
    local delta_x = xi_new - DESLOCA.bar[Paint.bar_sel].xi
    DESLOCA.bar[Paint.bar_sel].xi = xi_new
    if check_box_L_fixed:GetValue() == true then -- L fixed
      DESLOCA.bar[Paint.bar_sel].xf = DESLOCA.bar[Paint.bar_sel].xf + delta_x
      text_ctrl_xf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xf))
    else -- L free
      local L_new, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
      local angle_rad, angle_deg = Get_angle(S, C)
      spin_ctrl_angle:SetValue(angle_deg)
      DESLOCA.bar[Paint.bar_sel].L = L_new
      text_ctrl_L:ChangeValue(tostring(L_new))
    end
    DESLOCA.update_all()
  end
end)
-- Bar edit "yi"
panel:Connect(IDs.text_ctrl_yi, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local yi_new = tonumber(text_ctrl_yi:GetValue()) -- accept only numbers
  if text_ctrl_yi:GetValue() == "" or text_ctrl_yi:GetValue() == "-" or text_ctrl_yi:GetValue() == "." or text_ctrl_yi:GetValue() == "-." then
    -- wait number or error
  elseif yi_new == nil then
    wx.wxMessageBox("O valor inserido para \"yi\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  elseif Is_point(DESLOCA.bar[Paint.bar_sel].xi, yi_new, DESLOCA.bar[Paint.bar_sel].xf, DESLOCA.bar[Paint.bar_sel].yf, true) then
    wx.wxMessageBox("Barra muito pequena ou pontual, insira novas coordenadas.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
    text_ctrl_yi:ChangeValue("")
  else
    local delta_y = yi_new - DESLOCA.bar[Paint.bar_sel].yi
    DESLOCA.bar[Paint.bar_sel].yi = yi_new
    if check_box_L_fixed:GetValue() == true then -- L fixed
      DESLOCA.bar[Paint.bar_sel].yf = DESLOCA.bar[Paint.bar_sel].yf + delta_y
      text_ctrl_yf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yf))
    else -- L free
      local L_new, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
      local angle_rad, angle_deg = Get_angle(S, C)
      spin_ctrl_angle:SetValue(angle_deg)
      DESLOCA.bar[Paint.bar_sel].L = L_new
      text_ctrl_L:ChangeValue(tostring(L_new))
    end
    DESLOCA.update_all()
  end
end)

-- Bar edit "xf"
panel:Connect(IDs.text_ctrl_xf, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local xf_new = tonumber(text_ctrl_xf:GetValue()) -- accept only numbers
  if text_ctrl_xf:GetValue() == "" or text_ctrl_xf:GetValue() == "-" or text_ctrl_xf:GetValue() == "." or text_ctrl_xf:GetValue() == "-." then
    -- wait number or error
  elseif xf_new == nil then
    wx.wxMessageBox("O valor inserido para \"xf\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  elseif Is_point(DESLOCA.bar[Paint.bar_sel].xi, DESLOCA.bar[Paint.bar_sel].yi, xf_new, DESLOCA.bar[Paint.bar_sel].yf, true) then
    wx.wxMessageBox("Barra muito pequena ou pontual, insira novas coordenadas.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
    text_ctrl_xf:ChangeValue("")
  else
    local delta_x = xf_new - DESLOCA.bar[Paint.bar_sel].xf
    DESLOCA.bar[Paint.bar_sel].xf = xf_new
    if check_box_L_fixed:GetValue() == true then -- L fixed
      DESLOCA.bar[Paint.bar_sel].xi = DESLOCA.bar[Paint.bar_sel].xi + delta_x
      text_ctrl_xi:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xi))
    else -- L free
      local L_new, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
      local angle_rad, angle_deg = Get_angle(S, C)
      spin_ctrl_angle:SetValue(angle_deg)
      DESLOCA.bar[Paint.bar_sel].L = L_new
      text_ctrl_L:ChangeValue(tostring(L_new))
    end
    DESLOCA.update_all()
  end
end)
-- Bar edit "yf"
panel:Connect(IDs.text_ctrl_yf, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local yf_new = tonumber(text_ctrl_yf:GetValue()) -- accept only numbers
  if text_ctrl_yf:GetValue() == "" or text_ctrl_yf:GetValue() == "-" or text_ctrl_yf:GetValue() == "." or text_ctrl_yf:GetValue() == "-." then
    -- wait number or error
  elseif yf_new == nil then
    wx.wxMessageBox("O valor inserido para \"yf\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  elseif Is_point(DESLOCA.bar[Paint.bar_sel].xi, DESLOCA.bar[Paint.bar_sel].yi, DESLOCA.bar[Paint.bar_sel].xf, yf_new, true) then
    wx.wxMessageBox("Barra muito pequena ou pontual, insira novas coordenadas.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
    text_ctrl_yf:ChangeValue("")
  else
    local delta_y = yf_new - DESLOCA.bar[Paint.bar_sel].yf
    DESLOCA.bar[Paint.bar_sel].yf = yf_new
    if check_box_L_fixed:GetValue() == true then -- L fixed
      DESLOCA.bar[Paint.bar_sel].yi = DESLOCA.bar[Paint.bar_sel].yi + delta_y
      text_ctrl_yi:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yi))
    else -- L free
      local L_new, C, S = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
      local angle_rad, angle_deg = Get_angle(S, C)
      spin_ctrl_angle:SetValue(angle_deg)
      DESLOCA.bar[Paint.bar_sel].L = L_new
      text_ctrl_L:ChangeValue(tostring(L_new))
    end
    DESLOCA.update_all()
  end
end)

--[[ Bar edit "<" -- TODO: too many problems
panel:Connect(IDs.spin_ctrl_angle, wx.wxEVT_COMMAND_TEXT_UPDATED, function() -- this event happens too when wxEVT_SCROLL_LINEUP/DOWN happens
  local angle_new = tonumber(spin_ctrl_angle:GetValue()) -- accept only numbers
  if angle_new == nil or angle_new >= 361 or angle_new < -1 then -- sometimes it doesn't happen for some reason (using letters or big numbers)
    wx.wxMessageBox("O valor inserido para o ângulo \"<\" é inválido, insira apenas números entre 0 e 359.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    local angle_new_rad = (angle_new*2*math.pi)/360
    DESLOCA.bar[Paint.bar_sel].xf = DESLOCA.bar[Paint.bar_sel].xi + tonumber(text_ctrl_L:GetValue()) * math.cos(angle_new_rad)
    DESLOCA.bar[Paint.bar_sel].yf = DESLOCA.bar[Paint.bar_sel].yi + tonumber(text_ctrl_L:GetValue()) * math.sin(angle_new_rad)
    text_ctrl_xf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xf))
    text_ctrl_yf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yf))
    DESLOCA.eval_nodes()
    DESLOCA.all_matrix_update()
    DESLOCA.load_vector_update()
    DESLOCA.grid_fill_all()
    frame:Refresh()
  end
end)]]
-- Bar edit "<" (up)
panel:Connect(IDs.spin_ctrl_angle, wx.wxEVT_SCROLL_LINEUP, function()
  local angle_new = tonumber(spin_ctrl_angle:GetValue()) -- accept only numbers
  if angle_new == nil then
    wx.wxMessageBox("O valor inserido para o ângulo \"<\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    if angle_new < 359 then angle_new = angle_new + 1 else angle_new = 0 ; spin_ctrl_angle:SetValue(-1) end
    local angle_new_rad = (angle_new*2*math.pi)/360
    DESLOCA.bar[Paint.bar_sel].xf = DESLOCA.bar[Paint.bar_sel].xi + DESLOCA.bar[Paint.bar_sel].L * math.round(math.cos(angle_new_rad), 7)
    DESLOCA.bar[Paint.bar_sel].yf = DESLOCA.bar[Paint.bar_sel].yi + DESLOCA.bar[Paint.bar_sel].L * math.round(math.sin(angle_new_rad), 7)
    text_ctrl_xf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xf))
    text_ctrl_yf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yf))
    DESLOCA.update_all()
  end
end)
-- Bar edit "<" (down)
panel:Connect(IDs.spin_ctrl_angle, wx.wxEVT_SCROLL_LINEDOWN, function()
  local angle_new = tonumber(spin_ctrl_angle:GetValue()) -- accept only numbers
  if angle_new == nil then
    wx.wxMessageBox("O valor inserido para o Ã¢ngulo \"<\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    if angle_new > 0 then angle_new = angle_new - 1 else angle_new = 359 ; spin_ctrl_angle:SetValue(360) end
    local angle_new_rad = (angle_new*2*math.pi)/360
    DESLOCA.bar[Paint.bar_sel].xf = DESLOCA.bar[Paint.bar_sel].xi + DESLOCA.bar[Paint.bar_sel].L * math.round(math.cos(angle_new_rad), 7)
    DESLOCA.bar[Paint.bar_sel].yf = DESLOCA.bar[Paint.bar_sel].yi + DESLOCA.bar[Paint.bar_sel].L * math.round(math.sin(angle_new_rad), 7)
    text_ctrl_xf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].xf))
    text_ctrl_yf:ChangeValue(tostring(DESLOCA.bar[Paint.bar_sel].yf))
    DESLOCA.update_all()
  end
end)

-- Bar edit "E"
panel:Connect(IDs.text_ctrl_E, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local E_new = tonumber(text_ctrl_E:GetValue()) -- accept only numbers
  if text_ctrl_E:GetValue() == "" or text_ctrl_E:GetValue() == "." then
    -- wait number or error
  elseif E_new == nil then
    wx.wxMessageBox("O valor inserido para \"E\" é inválido, insira apenas números, e positivos.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].E = E_new
    DESLOCA.update_all()
  end
end)

-- Bar edit "A"
panel:Connect(IDs.text_ctrl_A, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local A_new = tonumber(text_ctrl_A:GetValue()) -- accept only numbers
  if text_ctrl_A:GetValue() == "" or text_ctrl_A:GetValue() == "." then
    -- wait number or error
  elseif A_new == nil then
    wx.wxMessageBox("O valor inserido para \"A\" é inválido, insira apenas números, e positivos.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else 
    DESLOCA.bar[Paint.bar_sel].A = A_new
    DESLOCA.update_all()
  end
end)

-- Bar edit "I"
panel:Connect(IDs.text_ctrl_I, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local I_new = tonumber(text_ctrl_I:GetValue()) -- accept only numbers
  if text_ctrl_I:GetValue() == "" or text_ctrl_I:GetValue() == "." then
    -- wait number or error
  elseif I_new == nil then
    wx.wxMessageBox("O valor inserido para \"I\" é inválido, insira apenas números, e positivos.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].I = I_new
    DESLOCA.update_all()
  end
end)

-- Bar distributed load "q"
panel:Connect(IDs.dist_load_text_ctrl_q, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local q_new = tonumber(dist_load_text_ctrl_q:GetValue()) -- accept only numbers
  if dist_load_text_ctrl_q:GetValue() == "" or dist_load_text_ctrl_q:GetValue() == "-" or dist_load_text_ctrl_q:GetValue() == "." or dist_load_text_ctrl_q:GetValue() == "-." then
    -- wait number or error
  elseif q_new == nil then
    wx.wxMessageBox("O valor inserido para \"q\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].dist.q = q_new
    DESLOCA.update_all()
  end
end)
-- Bar distributed load "local" option
panel:Connect(IDs.dist_load_option_local, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].dist.is_global = false
  DESLOCA.update_all()
end)
-- Bar distributed load "global em x" option
panel:Connect(IDs.dist_load_option_global_x, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].dist.is_global = true
  DESLOCA.bar[Paint.bar_sel].dist.global_dir = "x"
  DESLOCA.update_all()
end)
-- Bar distributed load "global em y" option
panel:Connect(IDs.dist_load_option_global_y, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].dist.is_global = true
  DESLOCA.bar[Paint.bar_sel].dist.global_dir = "y"
  DESLOCA.update_all()
end)

-- Bar nodal load "Fx" for i -- TODO: perhaps make a function that accept the type of force ("Fx", "Fy" or "Mz") and the node ("i" or "f") and automate this all saving some lines
panel:Connect(IDs.nodal_load_text_ctrl_Fx_i, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local Fxi_new = tonumber(nodal_load_text_ctrl_Fx_i:GetValue()) -- accept only numbers
  if nodal_load_text_ctrl_Fx_i:GetValue() == "" or nodal_load_text_ctrl_Fx_i:GetValue() == "-" or nodal_load_text_ctrl_Fx_i:GetValue() == "." or nodal_load_text_ctrl_Fx_i:GetValue() == "-." then
    -- wait number or error
  elseif Fxi_new == nil then
    wx.wxMessageBox("O valor inserido para \"Fx\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].i.Fx = Fxi_new
    DESLOCA.update_all()
  end
end)
-- Bar nodal load "Fy" for i
panel:Connect(IDs.nodal_load_text_ctrl_Fy_i, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local Fyi_new = tonumber(nodal_load_text_ctrl_Fy_i:GetValue()) -- accept only numbers
  if nodal_load_text_ctrl_Fy_i:GetValue() == "" or nodal_load_text_ctrl_Fy_i:GetValue() == "-" or nodal_load_text_ctrl_Fy_i:GetValue() == "." or nodal_load_text_ctrl_Fy_i:GetValue() == "-." then
    -- wait number or error
  elseif Fyi_new == nil then
    wx.wxMessageBox("O valor inserido para \"Fy\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].i.Fy = Fyi_new
    DESLOCA.update_all()
  end
end)
-- Bar nodal load "Mz" for i
panel:Connect(IDs.nodal_load_text_ctrl_Mz_i, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local Mzi_new = tonumber(nodal_load_text_ctrl_Mz_i:GetValue()) -- accept only numbers
  if nodal_load_text_ctrl_Mz_i:GetValue() == "" or nodal_load_text_ctrl_Mz_i:GetValue() == "-" or nodal_load_text_ctrl_Mz_i:GetValue() == "." or nodal_load_text_ctrl_Mz_i:GetValue() == "-." then
    -- wait number or error
  elseif Mzi_new == nil then
    wx.wxMessageBox("O valor inserido para \"Mz\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].i.Mz = Mzi_new
    DESLOCA.update_all()
  end
end)
-- Bar nodal load "local" option for i
panel:Connect(IDs.nodal_load_option_local_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.is_global = false
  DESLOCA.update_all()
end)
-- Bar nodal load "global" option for i
panel:Connect(IDs.nodal_load_option_global_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.is_global = true
  DESLOCA.update_all()
end)

-- Bar nodal load "Fx" for f -- TODO: perhaps make a function that accept the type of force ("Fx", "Fy" or "Mz") and the node ("i" or "f") and automate this all saving some lines
panel:Connect(IDs.nodal_load_text_ctrl_Fx_f, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local Fxf_new = tonumber(nodal_load_text_ctrl_Fx_f:GetValue()) -- accept only numbers
  if nodal_load_text_ctrl_Fx_f:GetValue() == "" or nodal_load_text_ctrl_Fx_f:GetValue() == "-" or nodal_load_text_ctrl_Fx_f:GetValue() == "." or nodal_load_text_ctrl_Fx_f:GetValue() == "-." then
    -- wait number or error
  elseif Fxf_new == nil then
    wx.wxMessageBox("O valor inserido para \"Fx\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].f.Fx = Fxf_new
    DESLOCA.update_all()
  end
end)
-- Bar nodal load "Fy" for f
panel:Connect(IDs.nodal_load_text_ctrl_Fy_f, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local Fyf_new = tonumber(nodal_load_text_ctrl_Fy_f:GetValue()) -- accept only numbers
  if nodal_load_text_ctrl_Fy_f:GetValue() == "" or nodal_load_text_ctrl_Fy_f:GetValue() == "-" or nodal_load_text_ctrl_Fy_f:GetValue() == "." or nodal_load_text_ctrl_Fy_f:GetValue() == "-." then
    -- wait number or error
  elseif Fyf_new == nil then
    wx.wxMessageBox("O valor inserido para \"Fy\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].f.Fy = Fyf_new
    DESLOCA.update_all()
  end
end)
-- Bar nodal load "Mz" for f
panel:Connect(IDs.nodal_load_text_ctrl_Mz_f, wx.wxEVT_COMMAND_TEXT_UPDATED, function()
  local Mzf_new = tonumber(nodal_load_text_ctrl_Mz_f:GetValue()) -- accept only numbers
  if nodal_load_text_ctrl_Mz_f:GetValue() == "" or nodal_load_text_ctrl_Mz_f:GetValue() == "-" or nodal_load_text_ctrl_Mz_f:GetValue() == "." or nodal_load_text_ctrl_Mz_f:GetValue() == "-." then
    -- wait number or error
  elseif Mzf_new == nil then
    wx.wxMessageBox("O valor inserido para \"Mz\" é inválido, insira apenas números.", "Erro!", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, panel)
  else
    DESLOCA.bar[Paint.bar_sel].f.Mz = Mzf_new
    DESLOCA.update_all()
  end
end)
-- Bar nodal load "local" option for f
panel:Connect(IDs.nodal_load_option_local_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.is_global = false
  DESLOCA.update_all()
end)
-- Bar nodal load "global" option for f
panel:Connect(IDs.nodal_load_option_global_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.is_global = true
  DESLOCA.update_all()
end)

-- Bar restrain "X" for i
panel:Connect(IDs.restrain_option_free_x_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.free_x = true
  DESLOCA.update_all()
end)
panel:Connect(IDs.restrain_option_fixed_x_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.free_x = false
  DESLOCA.update_all()
end)
-- Bar restrain "Y" for i
panel:Connect(IDs.restrain_option_free_y_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.free_y = true
  DESLOCA.update_all()
end)
panel:Connect(IDs.restrain_option_fixed_y_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.free_y = false
  DESLOCA.update_all()
end)
-- Bar restrain "Z" for i
panel:Connect(IDs.restrain_option_free_z_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.free_z = true
  DESLOCA.update_all()
end)
panel:Connect(IDs.restrain_option_fixed_z_i, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].i.free_z = false
  DESLOCA.update_all()
end)

-- Bar restrain "X" for f
panel:Connect(IDs.restrain_option_free_x_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.free_x = true
  DESLOCA.update_all()
end)
panel:Connect(IDs.restrain_option_fixed_x_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.free_x = false
  DESLOCA.update_all()
end)
-- Bar restrain "Y" for f
panel:Connect(IDs.restrain_option_free_y_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.free_y = true
  DESLOCA.update_all()
end)
panel:Connect(IDs.restrain_option_fixed_y_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.free_y = false
  DESLOCA.update_all()
end)
-- Bar restrain "Z" for f
panel:Connect(IDs.restrain_option_free_z_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.free_z = true
  DESLOCA.update_all()
end)
panel:Connect(IDs.restrain_option_fixed_z_f, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar[Paint.bar_sel].f.free_z = false
  DESLOCA.update_all()
end)

-- Bar display mode choice button
local bar_mode_choice_prev = 0
panel:Connect(IDs.bar_mode_choice, wx.wxEVT_COMMAND_CHOICE_SELECTED, function()
  local new_value = bar_mode_choice:GetCurrentSelection()
  if new_value == 0 and new_value ~= bar_mode_choice_prev then -- "Matrizes de rigidez e vetores de cargas nodais equivalentes" ; and to avoid calling the function if selected same option as current 
    DESLOCA.stiff_matrix_mode()
    DESLOCA.grid_fill_all()
    frame:Refresh()
  elseif new_value == 1 and new_value ~= bar_mode_choice_prev then -- "Forças nas extremidades em coordenadas locais" ; and to avoid calling the function if selected same option as current
    DESLOCA.forces_on_ends_mode()
    DESLOCA.grid_fill_all()
    frame:Refresh()
  end
  bar_mode_choice_prev = new_value
  DESLOCA.clear_status_bar(true, true)
end)

-- Bar rotation matrix option
panel:Connect(IDs.grid_rotation_mat_option_normal, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar_rotation_grid_fill()
  DESLOCA.clear_status_bar(true, true)
end)
panel:Connect(IDs.grid_rotation_mat_option_transposed, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.bar_rotation_grid_fill()
  DESLOCA.clear_status_bar(true, true)
end)

-- Structure display mode choice button
local structure_mode_choice_prev = 0
panel:Connect(IDs.structure_mode_choice, wx.wxEVT_COMMAND_CHOICE_SELECTED, function()
  local new_value = structure_mode_choice:GetCurrentSelection()
  if new_value ~= structure_mode_choice_prev then -- to avoid calling the function if selected same option as current
    DESLOCA.structure_mode_controls()
    frame:Refresh()
  end
  structure_mode_choice_prev = new_value
  DESLOCA.clear_status_bar(true, true)
end)

-- Deformed scale slider
panel:Connect(IDs.deformed_scale_slider, wx.wxEVT_SCROLL_THUMBTRACK, function()
  frame:Refresh()
  DESLOCA.clear_status_bar(true, true)
end)

-- Diagrams options
panel:Connect(IDs.diagram_option_N, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  frame:Refresh()
  DESLOCA.clear_status_bar(true, true)
end)
panel:Connect(IDs.diagram_option_Q, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  frame:Refresh()
  DESLOCA.clear_status_bar(true, true)
end)
panel:Connect(IDs.diagram_option_M, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  frame:Refresh()
  DESLOCA.clear_status_bar(true, true)
end)

-- Structure simplified option
panel:Connect(IDs.mat_estrut_simp_checkbox, wx.wxEVT_COMMAND_CHECKBOX_CLICKED, function()
  DESLOCA.grid_fill_all()
  frame:Refresh()
  DESLOCA.clear_status_bar(true, true)
end)

-- Structure load vector options
panel:Connect(IDs.grid_load_vector_structure_option1, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.structure_load_vector_grid_fill()
  DESLOCA.clear_status_bar(true, true)
end)
panel:Connect(IDs.grid_load_vector_structure_option2, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.structure_load_vector_grid_fill()
  DESLOCA.clear_status_bar(true, true)
end)
panel:Connect(IDs.grid_load_vector_structure_option3, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.structure_load_vector_grid_fill()
  DESLOCA.clear_status_bar(true, true)
end)

-- Results option
panel:Connect(IDs.grid_results_option_disp, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.results_grid_fill()
  DESLOCA.clear_status_bar(true, true)
end)
panel:Connect(IDs.grid_results_option_react, wx.wxEVT_COMMAND_RADIOBUTTON_SELECTED, function()
  DESLOCA.results_grid_fill()
  DESLOCA.clear_status_bar(true, true)
end)

-- Grid clicks
panel:Connect(IDs.grid_mat_local, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_mat_local) end)
panel:Connect(IDs.grid_load_vector_local, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_load_vector_local) end)
panel:Connect(IDs.grid_rotation_mat, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_rotation_mat) end)
panel:Connect(IDs.grid_mat_global, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_mat_global) end)
panel:Connect(IDs.grid_load_vector_global, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_load_vector_global) end)
panel:Connect(IDs.grid_f_ep, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_f_ep) end)
panel:Connect(IDs.grid_u, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_u) end)
panel:Connect(IDs.grid_f_local, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_f_local) end)
panel:Connect(IDs.grid_mat_estrut, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_mat_estrut) end)
panel:Connect(IDs.grid_load_vector_structure, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_load_vector_structure) end)
panel:Connect(IDs.grid_results, wx.wxEVT_GRID_CELL_LEFT_CLICK, function(event) DESLOCA.grid_click_handle(event, grid_results) end)

--- MENU EVENTS

-- Arquivo > Novo projeto
frame:Connect(IDs.file_menu_new_project, wx.wxEVT_COMMAND_MENU_SELECTED, function (event) -- TODO: make a window to warn the user the current structure will be deleted
  DESLOCA.bar = {
    [1] = {xi = 0, yi = 0, xf = 200, yf = 0, E = 20000, A = 100, I = 1000, global_node = {}, dist = {q = 0, is_global = false, global_dir = "x"},
           i = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = true, free_y = true, free_z = true},
           f = {Fx = 0, Fy = 0, Mz = 0, is_global = true, free_x = true, free_y = true, free_z = true}},
  }
  Paint.bar_sel = 1
  local L = DESLOCA.dimension(DESLOCA.bar[Paint.bar_sel])
  DESLOCA.bar[Paint.bar_sel].L = L -- to make sure
  DESLOCA.bar_choice_update()
  DESLOCA.update_all()
  DESLOCA.update_bar_info()
  frame:SetStatusText("Novo projeto iniciado!")
end)

-- Arquivo > Sair
frame:Connect(wx.wxID_EXIT, wx.wxEVT_COMMAND_MENU_SELECTED, function (event)
  frame:Close()
end)

-- Ajuda > Sobre
frame:Connect(wx.wxID_ABOUT, wx.wxEVT_COMMAND_MENU_SELECTED, function (event)
  wx.wxMessageBox("Este é o DESLOCA v"..DESLOCA.version..", uma ferramenta gráfico-didático-interativa para o ensino do Método dos Deslocamentos em estruturas bidimensionais. Está sendo desenvolvido em Lua, utilizando a biblioteca " .. wxlua.wxLUA_VERSION_STRING .. " para a construção da interface gráfica.",
                  "Sobre o DESLOCA",
                  wx.wxOK + wx.wxICON_INFORMATION,
                  frame)
end)

--############################################################################################################################################################
-- MAIN
--############################################################################################################################################################

-- Function to create a .ico file and set it as the program icon
function DESLOCA.icon()
  local file_name = "DESLOCA_icon.ico"
  
  local icon_file = assert(io.open(file_name, "wb"))
  local data = DESLOCA.icon_hex
  for j = 1, string.len(data)/2 do
    local char = string.char(tonumber(data:sub(2*j-1, 2*j), 16))
    icon_file:write(char) 
  end
  assert(icon_file:close())
  
  local icon = wx.wxIcon("DESLOCA_icon.ico", wx.wxBITMAP_TYPE_ICO)
  frame:SetIcon(icon)
  
  os.remove(file_name)
end

-- Init DESLOCA
function DESLOCA.init()
  
  -- Init the stiffness matrixes for the bar, since it creates global functions used later
  DESLOCA.stiff_matrix_mode()
  
  -- Init the whole core
  DESLOCA.update_all()
  
  -- Open DESLOCA window
  frame:Show(true)
  
  -- Set the program icon
  DESLOCA.icon()
  
  -- Small screen warning
  if SYSTEM.screen_is_small then
    wx.wxMessageBox(fmt("O DESLOCA é adequado para monitores com resolução %d x %d ou superior!\nO programa funcionará normalmente, mas talvez não seja possível observar todas as informações.", DESLOCA.width, DESLOCA.height), "Aviso", wx.wxOK + wx.wxICON_EXCLAMATION + wx.wxCENTRE, frame)
  end
  
  -- Set the welcome message in status bar
  status_bar:SetStatusText("Bem-vindo(a) ao DESLOCA v"..DESLOCA.version, 0)
end

DESLOCA.init()

-- Call wx.wxGetApp():MainLoop() last to start the wxWidgets event loop, otherwise the wxLua program will exit immediately
-- Does nothing if running from wxLua, wxLuaFreeze, or wxLuaEdit since the MainLoop is already running or will be started by the C++ program
wx.wxGetApp():MainLoop()

