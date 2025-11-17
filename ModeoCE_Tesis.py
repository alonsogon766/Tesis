# -*- coding: utf-8 -*-
"""
Created on Sat Sep 20 16:20:36 2025

@author: AlonsoGonzalez
"""
import pandas as pd
import os
from datetime import datetime
from gurobipy import Model, GRB, quicksum
import numpy as np

#Parámetros para escenarios
sv = 1  # Parámetro binario que habilita la generación solar en el vehículo. (1 = vehículo solar, 0 = vehículo eléctrico)
pf = 0  # Parámetro binario que determina el perfil fotovoltaico. (0 = calama, 1 = santiago, 2 = puerto montt)
year = 7 #Años 
usuarios = 4 #Cantidad de usuarios en la comunidad energética

ruta_salida = r"C:\Users\Administrator\Desktop\Alonso\Tesis\resultsCE\RSCE_CL_SV_7Y_4U_0.xlsx"

#Parámetros sujetos al periodo
if pf == 0:
    PVA = pd.read_excel('data/pv_profile_CALAMA.xlsx', header=None).iloc[:,0].to_numpy()         # Radiación incidente en el periodo t (Adimensional)  Data extraida de: https://nsrdb.nrel.gov/data-viewer
    PVA_sv = pd.read_excel('data/pv_profile_sv_CALAMA.xlsx', header=None).iloc[:,0].to_numpy()   # Radiación incidente en el vehículo solar en el periodo t (Adimensional)  
if pf == 1:
    PVA = pd.read_excel('data/pv_profile_SANTIAGO.xlsx', header=None).iloc[:,0].to_numpy()         # Radiación incidente en el periodo t (Adimensional)  Data extraida de: https://nsrdb.nrel.gov/data-viewer
    PVA_sv = pd.read_excel('data/pv_profile_sv_SANTIAGO.xlsx', header=None).iloc[:,0].to_numpy()   # Radiación incidente en el vehículo solar en el periodo t (Adimensional)  
if pf == 2:
    PVA = pd.read_excel('data/pv_profile_PUERTOMONTT.xlsx', header=None).iloc[:,0].to_numpy()         # Radiación incidente en el periodo t (Adimensional)  Data extraida de: https://nsrdb.nrel.gov/data-viewer
    PVA_sv = pd.read_excel('data/pv_profile_sv_PUERTOMONTT.xlsx', header=None).iloc[:,0].to_numpy()   # Radiación incidente en el vehículo solar en el periodo t (Adimensional)      
D   = pd.read_excel('dataCE/demandCE.xlsx', header=None).iloc[:, :usuarios].to_numpy().T              # Demanda del usuario i en la hora t (kWh) **MODIFICAR BASE DE DATOS**
Dsv  = pd.read_excel('dataCE/demand_svCE.xlsx',header=None).iloc[:, :usuarios].to_numpy().T                      # Demanda del vehículo solar/electrico en la hora t del usuario i  (km) **MODIFICAR BASE DE DATOS** Source: https://findingspress.org/article/10110
#CC_t = pd.read_excel('dataCE/price_buy.xlsx', header=None).iloc[:,0].to_numpy()                    # Precio de compra en el periodo t (u.m./kWh) **MODIFICAR BASE DE DATOS**
#CV_t = pd.read_excel('dataCE/price_sell.xlsx', header=None).iloc[:,0].to_numpy()                   # Precio de venta en el periodo t (u.m./kWh) **MODIFICAR BASE DE DATOS**
SVGRID = pd.read_excel('dataCE/sv_behaviorCE.xlsx', header=None).iloc[:, :usuarios].to_numpy().T                # Parámetro binario del comportamiento del vehículo solar (1: si está conectado, 0 e.o.c.)  **MODIFICAR BASE DE DATOS**
epsilon = pd.read_excel('data/epsilon.xlsx', header=None).iloc[:,0].to_numpy()                    # Porcentaje de degradación horaria de la batería, data extraida de paper de Abhi, suponiendo una degradación anual de un 1,2%
#Iveh= pd.read_excel('dataCE/Iveh.xlsx', header=None).iloc[usuarios:,0].to_numpy()                        # Parámetro binario si el usuario i tiene vehículo {1 si, 0 no}
#Isv = pd.read_excel('dataCE/Isv.xlsx', header=None).iloc[usuarios:,0].to_numpy()                        # Parámetro binario si el usuario i tiene vehículo solar, dado que tiene un vehículo {1 sv, 0 ev}

horas_a_simular = 24 * 365  # Simular una semana (168 horas)
PVA = PVA[:horas_a_simular]
D = D[:, :horas_a_simular]
Dsv = Dsv[:, :horas_a_simular]
SVGRID = SVGRID[:, :horas_a_simular]
epsilon = epsilon[:horas_a_simular]
epsilon = np.power(epsilon, year)
# ----------------------------------------------------

T = horas_a_simular 
I = usuarios


# Parámetros fijos
if sv == 1:
    DP = 0.7 #Depreciación horaria [u.m.]teniendo en cuenta el precio de $43500 y amortización de 7 años según el SII; Sources: https://aptera.us/reserve, https://www.sii.cl/valores_y_fechas/tabla_vida_util_activo_inmovilizado.html
if sv == 0:
    DP = 0.5 #Depreciación horaria [u.m.] teniendo en cuenta el precio de $30910 y amortización de 7 años según el SII; Sources: https://ev-database.org/car/2133/Renault-5-E-Tech-40kWh-95hp https://www.sii.cl/valores_y_fechas/tabla_vida_util_activo_inmovilizado.html
IPV    = 757.04   # Costo inversión PV (u.m./kW) Panel + inversor ##373.64 (u.m./kW) (panel) + 383.4 (u.m./kW) (inversor) = 757.04   https://www.solarstore.cl/producto/panel-solar-460w-48v-mono-dah-solar-hv ; https://www.solarstore.cl/producto/inversor-solis-s5-10kw-on-grid-hibrido-trifasico-certificado-sec
IBESS  = 384.75 # Costo inversión BESS (u.m./kWh) ##384.75                                                   https://www.solarstore.cl/producto/bateria-litio-pylontech-48kwh-48v-100ah-us-5000
PVsv = 0.7 #Potencia de los paneles fotovoltaicos del vehículo solar (kW)                                                   https://aptera.us/wp-content/uploads/2021/10/Specs-2022.pdf
KSV = 40 #Capacidad de la batería del vehículo solar. (kWh)                                                                  https://aptera.us/article/what-batteries-are-inside-aptera
CC= 0.24  #Precio fijado de compra de energía a la red (u.m./kWh)  **MODIFICAR**
CV= 0.11   #Precio fijado de venta de energía a la red (u.m./kWh)**MODIFICAR**
UB = 7.5   #Upper bound de potencia transmitida entre la red y el usuario. (kWh) o Potencia Instalada

CU = 20

#Parámetros de eficiencia
etacBESS = 0.9    #Eficiencia de carga de la BESS, data extraida de paper de Abhi
etadBESS = 0.9    #Eficiencia de descarga de la BESS, data extraida de paper de Abhi
etacSV   = 0.9    #Eficiencia de carga de la batería del SV, data extraida de paper de Abhi
etadSV   = 0.9    #Eficiencia de descarga de la batería del SV, data extraida de paper de Abhi
gamma_PV  = 1       #Factor de pérdida de eficiencia por inclinación para el panel fotovoltaico. (1: inclinación óptima)

if pf == 0:
    gamma_SV  = 0.981    #Factor de pérdida de eficiencia por inclinación para el panel fotovoltaico del SV.                    (Tener en cuenta que variará según la ubicación, para Calama= 0.981, Santiago = 0.934, Puerto Montt= 0.881)
if pf == 1:
    gamma_SV  = 0.934
if pf == 2:
    gamma_SV  = 0.881

if sv == 1:
    sigma = 0.075   #eficiencia del auto aptera 100Wh por milla, o 75Wh/km o 0.075kWh/km                                     #Source: https://www.youtube.com/watch?v=Kvpn7cteW5U
if sv == 0:
    sigma = 0.154   #eficiencia del auto renault 154Wh/km o 0.154kWh/km                                                      #Source: https://ev-database.org/car/2133/Renault-5-E-Tech-40kWh-95hp

#Parámetros degradación
alphaDDBESS = 0.000075   #Factor de degradación por deep discharge para la BESS (kWh)                                       , data extraida de paper de Abhi
alphaDDSV   = 0.000075   #Factor de degradación por deep discharge para la batería del vehículo solar (kWh)                    , data extraida de paper de Abhi
phi = 0.00022

##Deep dischargue establecido como máximo un 20% del total
betaSV   = 6   #capacidad a reducir por deep discharge para la batería del vehículo solar (kWh) como un 15% de la KSV #6

#Ramp up and down
CRate = 0.5  # Ratio de carga y descarga de la batería.

##BIG M
M = 10000 #Numero grande

f = 0.2  #Porcentaje de la batería mínimo 
h = 0.15  #Porcemtaje de la bateríaa reducir #0.15
r = 0.5  #Porcentaje de batería inicial

Delta = 1 #Tiempo de resolucion (h)

# Crear modelo Gurobi
m = Model('Modelo_TESIS')

##########################################################VARIABLES####################################################

# Variables de capacidad
pv = m.addVar(lb=0.0, name='pv')         # Capacidad PV instalada (kWp)
kBESS = m.addVar(lb=0.0, name='kBESS')   # Capacidad BESS instalada (kWh)

# Variables de flujo por periodo
rh = m.addVars(T, lb=0.0, name='rh_t')   # Energía transferida desde la red al hub de la comunidad energética. (kWh)
hr = m.addVars(T, lb=0.0, name='hr_t')   # Energía transferida desde desde el hub de la comunidad energética a la red. (kWh)
ph = m.addVars(T, lb=0.0, name='ph_t')   # Energía transferida desde la instalación fotovoltaica al hub de la comunidad energética. (kWh)
bh = m.addVars(T, lb=0.0, name='bh_t')   # Energía transferida desde la BESS al hub de la comunidad energética. (kWh)
hb = m.addVars(T, lb=0.0, name='hb_t')   # Energía transferida desde e. (kWh)
hd = m.addVars(I, T, lb=0.0, name='hd_it')   # Energía transferida desde el vehículo solar a la BESS en el periodo t. (kWh)
dh = m.addVars(I, T, lb=0.0, name='dh_it')   # Energía transferida desde la BESS al vehículo solar en el periodo t. (kWh)
ds = m.addVars(I, T, lb=0.0, name='ds_it')   # Energía transferida desde la BESS al usuario en el periodo t. (kWh)
sd = m.addVars(I, T, lb=0.0, name='sd_it')   # Energía transferida desde el sistema fotovoltaico al usuario en el periodo t. (kWh)

curtailment = m.addVars(T, lb=0.0, name= 'curt_t')

# Variables de estado de carga
socBESS = m.addVars(T, lb=0.0, name='socBESS_t') #Estado de carga de la BESS en el periodo t. (kWh)
socSV   = m.addVars(I, T, lb=0.0, name='socSV_it')   #Estado de carga de la batería del vehículo solar del usuario i en el periodo t. (kWh)

# Binarias
yBESS = m.addVars(T, vtype=GRB.BINARY, name='yBESS_t')      # {1, si la BESS realiza un deep discharge, 0 e.o.c.}
ySV   = m.addVars(I,T, vtype=GRB.BINARY, name='ySV_t')        # {1, si la batería del vehículo solar realiza un deep discharge, 0 e.o.c.}

#Variables de degradación Demanda no servida
betaBESS = m.addVars(T, lb=0.0, name='betaBESS_t')   #capacidad a reducir por deep discharge para la BESS en el periodo t                               establecido en la restricción 12

zBESS = m.addVars(T, vtype=GRB.BINARY, name='zBESS_t')    # {1 si BESS esta en modo CARGA, 0 si esta en modo DESCARGA}
zSV   = m.addVars(I, T, vtype=GRB.BINARY, name='zSV_it')   # {1 si SV[i] esta en modo CARGA, 0 si esta en modo DESCARGA}

#Degradación por ciclos
cum_flow_BESS = m.addVars(T, lb=0.0, name='cum_flow_BESS')
cum_flow_SV   = m.addVars(I, T, lb=0.0, name='cum_flow_SV')
n_cycles_SV   = m.addVars(I, T, vtype=GRB.INTEGER, lb=0.0, name='n_cycles_SV')

#Auxiliares
sum_yBESS = m.addVars(T, lb=0.0, name='sum_yBESS_t')
sum_ySV = m.addVars(I, T, lb=0.0, name='sum_ySV_it')
############################################################FUNCIÓN OBJETIVO############################################
Investment_Cost = pv * (IPV) + kBESS * (IBESS)
Annual_Grid_Cost = quicksum(CC * rh[t] - CV * hr[t] + curtailment[t] * CU for t in range(T))
#Annual_Grid_Cost = quicksum(CC * rh[t] - CV * hr[t] for t in range(T))
Annual_Depreciation_Cost = DP * T * usuarios
depreciation_years = min(year, 7)

m.setObjective(Investment_Cost + year * Annual_Grid_Cost + depreciation_years * Annual_Depreciation_Cost, GRB.MINIMIZE)
###########################################################RESTRICCIONES###############################################

# (2) Disponibilidad fotovoltaica
# ph[t] == PVAt[t] * gamma_PV * pv * PVF    for all t
for t in range(T):
    m.addConstr( ph[t] + curtailment[t] == PVA[t] * gamma_PV * pv * Delta, name=f'cons2_t{t}')
    #m.addConstr( ph[t] == PVA[t] * gamma_PV * pv * Delta, name=f'cons2_t{t}')

# (3) Estado de la carga de la BESS
# ηc_BESS * hb[t] + socBESS[t-1] == (1/ηd_BESS) * bh[t] + socBESS[t]  for t>=1
for t in range(1, T):
    lhs = etacBESS * hb[t] + socBESS[t-1]
    rhs = (1.0/etadBESS) * bh[t] + socBESS[t]
    m.addConstr( lhs == rhs, name=f'cons3_t{t}')

# (4) Condición de borde de la BESS 
m.addConstr( socBESS[0] == kBESS*r, name='cons4')

# (5) Balance HUBVG
# rh[t] + ph[t] + bh[t] + quicksum(dh[i,t]  for i in range(I)) = hr[t] + hb[t] + quicksum (hd[i,t] for i in range(I))
for t in range(T):
    m.addConstr( rh[t] + ph[t] + bh[t] + quicksum(dh[i,t]  for i in range(I)) == hr[t] + hb[t] + quicksum (hd[i,t] for i in range(I)), name=f'cons6_t{t}')

# (6) Demanda usuario 
# hd[i,t] + sd[i,t] * SVGRID[i,t] == D[i,t] + dh[i,t] + ds[i,t] * SVGRID[i,t] 
for t in range(T):
    for i in range(I):
        m.addConstr( hd[i,t] + sd[i,t] * SVGRID[i,t] == D[i,t] + dh[i,t] + ds[i,t] * SVGRID[i,t], name=f'cons7_i{i}_t{t}')

m.addConstr(cum_flow_BESS[0] == etacBESS * hb[0], name='cum_bess_0')
for t in range(1, T):
    m.addConstr(cum_flow_BESS[t] == cum_flow_BESS[t-1] + etacBESS * hb[t], name=f'cum_bess_t{t}')
    
# (7) Estado de carga y capacidad de la BESS
# soct_BESS <= kBESS * (epsilon[t]) - alpha_DD * sum_{t'<=t} y_t' - cumflow_BESS*phi   for all t
for t in range(T):
    cap_time = kBESS * epsilon[t]
    depdd = alphaDDBESS * sum_yBESS[t]
    deg_ciclo = cum_flow_BESS[t] * phi
    
    m.addConstr(socBESS[t] <= cap_time - depdd - deg_ciclo, name=f'cons7_cap_BESS_t{t}')
    m.addConstr( socBESS[t] <= (kBESS * (epsilon[t])) - alphaDDBESS * sum_yBESS[t], name=f'cons7_t{t}')
# (NUEVO) Definición de la suma acumulativa de 'yBESS'
m.addConstr(sum_yBESS[0] == yBESS[0], name='cons_sum_yBESS_0')
for t in range(1, T):
        m.addConstr(sum_yBESS[t] == sum_yBESS[t-1] + yBESS[t], name=f'cons_sum_yBESS_t{t}')

# (8) Estado de carga del vehículo solar 
# ηc_SV * SVGRID[i,t] * ds[i,t] + socSV[t-1] + PVsv * gamma_SV * PVA[t] == (1/ηd_SV) * SVGRID[it] * sd[i,t] + Dsv[i,t] + socSV[t]   for t>=2
for t in range(1, T):
    #Con generación solar
    for i in range(I):
        if sv == 1:
            lhs = etacSV * SVGRID[i,t] * ds[i,t] + socSV[i,t-1] + PVsv * gamma_SV * PVA[t]
            #Sin generación solar
        if sv == 0:
            lhs = etacSV * SVGRID[i,t] * ds[i,t] + socSV[i,t-1]
        rhs = (1.0/etadSV) * SVGRID[i,t] * sd[i,t] + Dsv[i,t] * sigma + socSV[i,t]
        m.addConstr( lhs == rhs, name=f'cons8_i{i}_t{t}')
   
# (9) Condición de borde de la batería del vehículo solar
for i in range(I):
    m.addConstr(socSV[i,0] == KSV*r, name='cons9_socSV_i{i}')


for i in range(I):
    input_sv_0 = etacSV * SVGRID[i,0] * ds[i,0] 
    m.addConstr(cum_flow_SV[i,0] == input_sv_0, name=f'cum_sv_i{i}_0')
    for t in range(1, T):
        input_sv_t = etacSV * SVGRID[i,t] * ds[i,t]
        m.addConstr(cum_flow_SV[i,t] == cum_flow_SV[i,t-1] + input_sv_t, name=f'cum_sv_i{i}_t{t}')
epsilon_step = 0.0001
for t in range(T):
    for i in range(I):
        m.addConstr(KSV * n_cycles_SV[i,t] >= cum_flow_SV[i,t] - KSV + epsilon_step, name=f'step_trigger_sv_i{i}_t{t}') 
        
# (10) Estado de carga y capacidad máxima de la batería del vehículo solar 
# socSV[i,t] <= KSV * (epsilon[t]) - alphaDDSV * sum_{t'<=t} ySV[i,t']
for t in range(T):
    for i in range(I):
        cap_time_sv = KSV * epsilon[t]
        deg_abuso_sv = alphaDDSV * sum_ySV[i,t]
        deg_ciclo_sv = n_cycles_SV[i,t] * (KSV * phi)
        m.addConstr(socSV[i,t] <= cap_time_sv - deg_abuso_sv - deg_ciclo_sv, name=f'cons10_cap_SV_t{t}')
# (NUEVO) Definición de la suma acumulativa de 'ySV'
for i in range(I):
    m.addConstr(sum_ySV[i,0] == ySV[i,0], name=f'cons_sum_ySV_i{i}_0')
    for t in range(1, T):
        m.addConstr(sum_ySV[i,t] == sum_ySV[i,t-1] + ySV[i,t], name=f'cons_sum_ySV_i{i}_t{t}')
        
# (11) Cantidad mínima de la batería y deep discharge del vehículo solar
# KSV * f * (epsilon[t]) - betaSV * ySV[i,t] <= socSV[i,t]  for all t>1
for t in range(1, T):
    for i in range(I):
        m.addConstr(KSV * f * (epsilon[t])  - betaSV * ySV[i,t] <= socSV[i,t], name=f'cons11_low_i{i}_t{t}')

# (12) Límite Mínimo del Estado de Carga (SOC) de la BESS
# kBESS * f - betaBESS <= socBESS[t] for all t>1
for t in range(1, T):
    m.addConstr(socBESS[t] >= kBESS * f - betaBESS[t], name=f'cons12_low_linear_t{t}')
    
# (13) Límite superior de activación de la penalización.
# betaBESSt <= M * yBESS[t]  for all t>1
for t in range(1, T):
    m.addConstr(betaBESS[t] <= M * yBESS[t])
    
# (14) Magnitud máxima de la penalización.    
# betaBESS[t] <= kBESS * h 
for t in range(1, T):
    m.addConstr(betaBESS[t] <= kBESS * h)
    
# (15) Límite inferior de activación de la penalización.
# betaBESS[t] >= kBESS * h - M (1 - yBESS[t])
for t in range(1, T):
    m.addConstr(betaBESS[t] >= kBESS * h - M * (1 - yBESS[t]))
    
# (16) Limite de transferencia desde el usuario a la red  
# hr[t] <= UB for all t
for t in range(T):
    m.addConstr(hr[t] <= UB * usuarios, name = f'cons16_usgridUB{t}')

# (17) Limite de transferencia desde la red al usuario 
# rh[t] <= UB for all t 
for t in range(T):
    m.addConstr(rh[t] <= UB * usuarios, name = f'cons17_gridusUB{t}')
    
# (18) y (19) Limite de inyección y extracción de la BESS respectivamente
# socBESS[t] - socBESS[t-1] <= kBESS * CRate && socBESS[t-1] - socBESS[t] <= kBESS * CRate
for t in range(1, T):
    m.addConstr(socBESS[t] - socBESS[t-1] <= kBESS * CRate, name=f'18ramp_up_soc_t{t}')
    m.addConstr(socBESS[t-1] - socBESS[t] <= kBESS * CRate, name=f'19ramp_down_soc_t{t}')
    
# (20) y (21) Limite de inyección y extracción de vehículo solar respectivamente
# socSV[i,t] - socSV[i,t-1] <= KSV * CRate && socSV[i,t-1] - socSV[i,t] <= KSV * CRate
for t in range(1, T):
    for i in range (I):
        m.addConstr(socSV[i,t] - socSV[i,t-1] <= KSV * CRate, name=f'20ramp_up_soc_t{t}')
        m.addConstr(socSV[i,t-1] - socSV[i,t] <= KSV * CRate, name=f'21ramp_down_soc_t{t}')    

# (22) Condición de Borde Final BESS 
# socBESS[T-1] >= socBESS[0] for t=|T|
m.addConstr(socBESS[T-1] >= socBESS[0], name='cons_final_soc_bess')

# (23) Condición de Borde Final SV
# soc_sv[T-1] >= soc_sv[0] for t=|T|
for i in range(I):
    m.addConstr(socSV[i,T-1] >= socSV[i,0], name='cons_final_soc_sv')

# (24) y (25) O Carga O Descarga BESS (Mutuamente Excluyente); zBESS = 1 -> Modo CARGA (hb > 0, bh = 0); zBESS = 0 -> Modo DESCARGA (hb = 0, bh > 0)
for t in range(T):
    # (24) Carga (hb) solo permitida si zBESS = 1
    m.addConstr(hb[t] <= M * zBESS[t], name=f'cons24_BESS_solo_carga_t{t}')
    
    # (25) Descarga (bh) solo permitida si zBESS = 0
    m.addConstr(bh[t] <= M * (1 - zBESS[t]), name=f'cons25_BESS_solo_descarga_t{t}')

# (26) y (27) O Carga O Descarga SV (Mutuamente Excluyente) ; zSV = 1 -> Modo CARGA (hd > 0 o ds > 0, dh = 0 y sd = 0) ; zSV = 0 -> Modo DESCARGA (hd = 0 y ds = 0, dh > 0 o sd > 0)
for t in range(T):
    for i in range(I):
        # (26) Flujos de CARGA (hd, ds) solo permitidos si zSV[i,t] = 1
        m.addConstr(hd[i,t] <= M * zSV[i,t], name=f'cons26_SV_solo_carga_G2V_i{i}_t{t}')
        m.addConstr(ds[i,t] <= M * zSV[i,t], name=f'cons26_SV_solo_carga_H2V_i{i}_t{t}')
        
        # (27) Flujos de DESCARGA (dh, sd) solo permitidos si zSV[i,t] = 0
        m.addConstr(dh[i,t] <= M * (1 - zSV[i,t]), name=f'cons27_SV_solo_descarga_V2G_i{i}_t{t}')
        m.addConstr(sd[i,t] <= M * (1 - zSV[i,t]), name=f'cons27_SV_solo_descarga_V2H_i{i}_t{t}')

# -------------------------
# Parámetros de solver y optimización
# -------------------------
m.Params.OutputFlag = 1    # 1 para ver salida, 0 para silencio
m.Params.DualReductions = 1
m.Params.MIPGap = 0.01
m.Params.Threads = 0
m.Params.Heuristics = 0.1
m.Params.NodefileStart = 100.0
m.Params.Method = 1
# -------------------------
# Optimizar
# -------------------------
m.optimize()

import pandas as pd

# Crear carpeta de salida si no existe
#ruta_salida = r"D:\Programacion\PYTHON\Tesis\resultsCE\RSCE_CL_NOSV_1Y_3U.xlsx"
os.makedirs(os.path.dirname(ruta_salida), exist_ok=True)

# ==========================================================
# VERIFICAR SOLUCION
# ==========================================================

if m.status == GRB.OPTIMAL:
    print("\n=== SOLUCION OPTIMA ENCONTRADA ===\n")
    print(f"Valor funcion objetivo: {m.objVal:.2f} [u.m.]")
    
    # =====================================================
    # PRE-CALCULO de demandas totales para graficos
    # =====================================================
    total_D_hogar_t = D.sum(axis=0) 
    total_D_sv_t = (Dsv * sigma).sum(axis=0) 
    
    # =====================================================
    # Variables comunes
    # =====================================================
    pv_opt = m.getVarByName("pv").X
    kBESS_opt = m.getVarByName("kBESS").X
    print(f"Capacidad PV instalada: {pv_opt:.2f} kW")
    print(f"Capacidad BESS instalada (Nominal): {kBESS_opt:.2f} kWh")

    # Calcular degradacion BESS
    total_deep_discharges_BESS = sum([m.getVarByName(f"yBESS_t[{t}]").X for t in range(T)])
    final_epsilon = epsilon[T-1]    
    
    capacidad_final_degradada_BESS = (kBESS_opt * final_epsilon) - (alphaDDBESS * total_deep_discharges_BESS)
    if capacidad_final_degradada_BESS < 0:
        capacidad_final_degradada_BESS = 0.0
    print(f"Capacidad BESS Final (Degradada): {capacidad_final_degradada_BESS:.2f} kWh")

    # =====================================================
    # NUEVO: Calculo de Indicadores Clave (KPIs)
    # =====================================================
    
    # --- KPIs Comunales (Anuales, basados en T=8760h) ---
    total_rh = sum(m.getVarByName(f"rh_t[{t}]").X for t in range(T)) # Compra Red
    total_hr = sum(m.getVarByName(f"hr_t[{t}]").X for t in range(T)) # Venta Red
    total_ph = sum(m.getVarByName(f"ph_t[{t}]").X for t in range(T)) # Gen PV
    total_hb = sum(m.getVarByName(f"hb_t[{t}]").X for t in range(T)) # Carga BESS
    total_bh = sum(m.getVarByName(f"bh_t[{t}]").X for t in range(T)) # Descarga BESS
    
    total_D_hogar_anual = D.sum() # Demanda Hogares (anual)
    total_D_sv_anual = (Dsv * sigma).sum() # Demanda Vehiculos (anual)
    total_consumo_comunidad = total_D_hogar_anual + total_D_sv_anual

    # Self-Sufficiency Ratio (SSR): % de demanda cubierto por generacion propia (PV)
    ssr = 0.0
    if total_consumo_comunidad > 0.01:
        ssr = (total_ph / total_consumo_comunidad) * 100.0
        
    # Self-Consumption Ratio (SCR): % de PV generada que se consume/almacena (no se vende)
    scr = 0.0
    if total_ph > 0.01:
        scr = ((total_ph - total_hr) / total_ph) * 100.0

    # Costos
    investment_cost_total = (pv_opt * IPV) + (kBESS_opt * IBESS)
    annual_grid_cost_neto = (CC * total_rh) - (CV * total_hr) # Costo neto de red (1 año)
    annual_depreciation_cost = DP * T * usuarios # Costo depreciacion (1 año)
    depreciation_years_applied = min(year, 7)
    total_depreciation_cost_project = annual_depreciation_cost * depreciation_years_applied
    
    # =====================================================
    # Data resumida (KPIs)
    # =====================================================
    datos_resumen = {
        'Metrica': [
            # --- Costos del Proyecto ---
            'Valor Funcion Objetivo (u.m.)',
            'Costo Inversion Inicial (u.m.)',
            'Costo Total Depreciacion Vehiculos (u.m.)',
            f'Años de Depreciacion Aplicados',
            f'Costo Neto Red (Total {year} años) (u.m.)',
            # --- KPIs Anuales Comunidad ---
            'Costo Neto Red (1 Año) (u.m.)',
            'Total Compra Red (kWh/año)',
            'Total Venta Red (kWh/año)',
            'Total Generacion PV (kWh/año)',
            'Total Consumo Hogares (kWh/año)',
            'Total Consumo Vehiculos (kWh/año)',
            'Ratio de Autosuficiencia (SSR) (%)',
            'Ratio de Autoconsumo (SCR) (%)',
            # --- Inversion Comunidad ---
            'Capacidad PV Instalada (kW)',
            'Capacidad BESS Nominal (kWh)',
            'Capacidad BESS Final (Degradada) (kWh)',
            'Total Deep Discharges BESS (ciclos/año)'
        ],
        'Valor': [
            # --- Costos del Proyecto ---
            m.objVal,
            investment_cost_total,
            total_depreciation_cost_project,
            depreciation_years_applied,
            annual_grid_cost_neto * year,
            # --- KPIs Anuales Comunidad ---
            annual_grid_cost_neto,
            total_rh,
            total_hr,
            total_ph,
            total_D_hogar_anual,
            total_D_sv_anual,
            ssr,
            scr,
            # --- Inversion Comunidad ---
            pv_opt,
            kBESS_opt,
            capacidad_final_degradada_BESS,
            total_deep_discharges_BESS
        ]
    }

    # =====================================================
    # Data SV (Se añade a datos_resumen)
    # =====================================================
    for i in range(I):
        total_dd_sv_i = sum([m.getVarByName(f"ySV_t[{i},{t}]").X for t in range(T)])
        
        cap_final_sv_i = (KSV * final_epsilon) - (alphaDDSV * total_dd_sv_i)
        if cap_final_sv_i < 0:
            cap_final_sv_i = 0.0
        
        print(f"Capacidad SV Final Usuario {i} (Degradada): {cap_final_sv_i:.2f} kWh")

        # KPIs por usuario
        user_total_D_hogar = D[i, :].sum()
        user_total_D_sv = (Dsv[i, :] * sigma).sum()
        user_total_sd_V2H = sum(m.getVarByName(f"sd_it[{i},{t}]").X for t in range(T)) # V2H (Vehiculo a Casa)
        user_total_ds_H2V = sum(m.getVarByName(f"ds_it[{i},{t}]").X for t in range(T)) # H2V (Casa a Vehiculo)
        user_total_dh_V2G = sum(m.getVarByName(f"dh_it[{i},{t}]").X for t in range(T)) # V2G (Vehiculo a Hub)
        user_total_hd_G2V = sum(m.getVarByName(f"hd_it[{i},{t}]").X for t in range(T)) # G2V (Hub a Vehiculo)
        user_total_charge_sv = user_total_ds_H2V + user_total_hd_G2V
        user_total_discharge_sv = user_total_sd_V2H + user_total_dh_V2G

        # Añadir separador y datos de este usuario
        datos_resumen['Metrica'].append(f'--- Metricas Usuario {i} ---')
        datos_resumen['Valor'].append(pd.NA) 
        
        datos_resumen['Metrica'].append(f'Usuario {i} - Consumo Hogar (kWh/año)')
        datos_resumen['Valor'].append(user_total_D_hogar)
        datos_resumen['Metrica'].append(f'Usuario {i} - Consumo Vehiculo (kWh/año)')
        datos_resumen['Valor'].append(user_total_D_sv)
        
        datos_resumen['Metrica'].append(f'Usuario {i} - Carga Total Vehiculo (kWh/año)')
        datos_resumen['Valor'].append(user_total_charge_sv)
        datos_resumen['Metrica'].append(f'Usuario {i} - (Detalle) Carga G2V (Hub) (kWh/año)')
        datos_resumen['Valor'].append(user_total_hd_G2V)
        
        datos_resumen['Metrica'].append(f'Usuario {i} - Descarga Total Vehiculo (kWh/año)')
        datos_resumen['Valor'].append(user_total_discharge_sv)
        datos_resumen['Metrica'].append(f'Usuario {i} - (Detalle) Descarga V2G (Hub) (kWh/año)')
        datos_resumen['Valor'].append(user_total_dh_V2G)
        datos_resumen['Metrica'].append(f'Usuario {i} - (Detalle) Descarga V2H (Casa) (kWh/año)')
        datos_resumen['Valor'].append(user_total_sd_V2H)

        datos_resumen['Metrica'].append(f'Usuario {i} - Capacidad SV Final (Degradada) (kWh)')
        datos_resumen['Valor'].append(cap_final_sv_i)
        
        datos_resumen['Metrica'].append(f'Usuario {i} - Total Deep Discharges SV (ciclos/año)')
        datos_resumen['Valor'].append(total_dd_sv_i)

    # !! CORREGIDO: El DataFrame se crea DESPUES de que el bucle for llene los datos
    df_resumen = pd.DataFrame(datos_resumen)

    
    # ======================================================
    # EXPORTAR RESULTADOS A EXCEL
    # ======================================================
    with pd.ExcelWriter(ruta_salida, engine='xlsxwriter') as writer:
        df_resumen.to_excel(writer,index=False,sheet_name='Resumen_Inversion')
        worksheet_resumen = writer.sheets['Resumen_Inversion']
        formato_numerico_resumen = writer.book.add_format({'num_format': '#,##0.00', 'bold': True})
        formato_texto_resumen = writer.book.add_format({'bold': False})
        formato_separador = writer.book.add_format({'bold': True, 'bg_color': '#D3D3D3'})
        
        # Aplicar formato condicional para separadores
        worksheet_resumen.conditional_format('A1:B100', {
            'type': 'text',
            'criteria': 'containing',
            'value': '---',
            'format': formato_separador
        })
        
        worksheet_resumen.set_column('A:A', 50, formato_texto_resumen) # Columna mas ancha
        worksheet_resumen.set_column('B:B', 20, formato_numerico_resumen)


        # --- Hojas por Usuario (Una por cada 'i') ---
        for i in range(I):
            sheet_name = f'Resultados_Usuario_{i}'
            print(f"Generando hoja: {sheet_name}...")

            # 1. DATOS COMUNALES
            resultados_i = pd.DataFrame({
                'Periodo': list(range(1, T + 1)),
                'PVA(kW/m2)': PVA,
                'rh_Comunal(kWh)': [m.getVarByName(f"rh_t[{t}]").X for t in range(T)],
                'hr_Comunal(kWh)': [m.getVarByName(f"hr_t[{t}]").X for t in range(T)],
                'ph_Comunal(kWh)': [m.getVarByName(f"ph_t[{t}]").X for t in range(T)],
                'bh_Comunal(kWh)': [m.getVarByName(f"bh_t[{t}]").X for t in range(T)],
                'hb_Comunal(kWh)': [m.getVarByName(f"hb_t[{t}]").X for t in range(T)],
                'socBESS(kWh)': [m.getVarByName(f"socBESS_t[{t}]").X for t in range(T)],
                'yBESS(bin)': [m.getVarByName(f"yBESS_t[{t}]").X for t in range(T)],
                'betaBESS(kWh)': [m.getVarByName(f"betaBESS_t[{t}]").X for t in range(T)],
                # NUEVO: Columnas de demanda total (para graficos comunales)
                'Total_Hogar_Demand(kWh)': total_D_hogar_t,
                'Total_SV_Demand(kWh)': total_D_sv_t,
            })

            # 2. DATOS POR USUARIO          
            # Manejar Dsv (que es 2D, forma [I, T])
            dsv_usuario = Dsv[i,:] if Dsv.ndim > 1 else Dsv
            resultados_i[f'D_sv_i{i}(kWh)'] = dsv_usuario * sigma # Consumo en kWh
            resultados_i[f'D_Hogar_i{i}(kWh)'] = D[i, :]

            resultados_i[f'hd_i{i}(kWh)'] = [m.getVarByName(f"hd_it[{i},{t}]").X for t in range(T)]
            resultados_i[f'dh_i{i}(kWh)'] = [m.getVarByName(f"dh_it[{i},{t}]").X for t in range(T)]
            resultados_i[f'ds_i{i}(kWh)'] = [m.getVarByName(f"ds_it[{i},{t}]").X for t in range(T)]
            resultados_i[f'sd_i{i}(kWh)'] = [m.getVarByName(f"sd_it[{i},{t}]").X for t in range(T)]
            resultados_i[f'socSV_i{i}(kWh)'] = [m.getVarByName(f"socSV_it[{i},{t}]").X for t in range(T)]
            resultados_i[f'ySV_i{i}(bin)'] = [m.getVarByName(f"ySV_t[{i},{t}]").X for t in range(T)]

            # -----------------------------------------------
            # Identificar nuevas columnas por indice (A=0)
            # -----------------------------------------------
            # A-J (10) son comunales
            # K: Total_Hogar_Demand (idx 10)
            # L: Total_SV_Demand (idx 11)
            # M: D_sv_i (idx 12)
            # N: D_Hogar_i (idx 13)
            # O: hd_i (idx 14)
            # P: dh_i (idx 15)
            # Q: ds_i (idx 16)
            # R: sd_i (idx 17)
            # S: socSV_i (idx 18)
            # T: ySV_i (idx 19)
            # Total 20 columnas (A-T)
            
            col_periodo = '$A'    # A
            col_rh_comunal = '$C' # C
            col_hr_comunal = '$D' # D
            col_ph_comunal = '$E' # E
            col_socBESS = '$H'    # H
            col_total_hogar_D = '$K' # K
            col_total_sv_D = '$L'    # L
            col_user_sv_D = '$M'     # M
            col_user_hogar_D = '$N'  # N
            col_user_socSV = '$S'    # S
            

            sum_row_data = {}
            for col in resultados_i.columns:
                if col == 'Periodo':
                    sum_row_data[col] = 'TOTAL'
                elif 'soc' in col.lower() or 'pva' in col.lower() or 'beta' in col.lower():
                    sum_row_data[col] = np.nan 
                else:
                    sum_row_data[col] = resultados_i[col].sum()
            sum_row_df = pd.DataFrame(sum_row_data, index=[T]) # Ponerlo al final
            resultados_con_total = pd.concat([resultados_i, sum_row_df])
            resultados_con_total.to_excel(writer, index=False, sheet_name=sheet_name)
            worksheet_horarios = writer.sheets[sheet_name]
            formato_numerico = writer.book.add_format({'num_format': '#,##0.000'})
            formato_total_fila = writer.book.add_format({
                'bold': True,
                'num_format': '#,##0.000',
                'top': 1 
            })

            # Formatear columnas (ahora A-T)
            for idx, col in enumerate(resultados_con_total.columns):
                if col == 'Periodo':
                    worksheet_horarios.set_column(idx, idx, 10)
                else:
                    worksheet_horarios.set_column(idx, idx, 18, formato_numerico)

            worksheet_horarios.set_row(T + 1, None, formato_total_fila)
            # Obtener el workbook
            workbook = writer.book

            data_rows = T + 1
            
            chart_x_scale = 1.8
            chart_y_scale = 1.3
            
            # =================================================
            # GRAFICO 1 Demanda hogar vs SOC vehiculo
            # =================================================
            chart1 = workbook.add_chart({'type': 'line'})

            # Añadir serie de Demanda Hogar (Col N)
            chart1.add_series({
                'name':f"='{sheet_name}'!{col_user_hogar_D}$1",#'D_Hogar_i{i}(kWh)'
                'categories':f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",# Eje X (Periodo)
                'values':f"='{sheet_name}'!{col_user_hogar_D}$2:{col_user_hogar_D}${data_rows}", # Eje Y (Demanda)
                'y2_axis': 0, # Eje Y primario
            })

            # Añadir serie de SOC del Vehiculo (Col S)
            chart1.add_series({
                'name': f"='{sheet_name}'!{col_user_socSV}$1", # 'socSV_i{i}(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}", # Eje X (Periodo)
                'values': f"='{sheet_name}'!{col_user_socSV}$2:{col_user_socSV}${data_rows}", # Eje Y (socSV)
                'y2_axis': 1, # Eje Y secundario
            })

            chart1.set_title({'name': f'User {i} Home Demand vs. Vehicle SOC'})
            chart1.set_x_axis({'name': 'Period (hour)'})
            chart1.set_y_axis({'name': 'Home Demand (kWh)'})
            chart1.set_y2_axis({'name': 'Vehicle SOC (kWh)'})

            worksheet_horarios.insert_chart('U2', chart1, {'x_scale': chart_x_scale, 'y_scale': chart_y_scale})

            # =================================================
            # GRAFICO 2 - Demandas Usuario
            # =================================================
            chart2 = workbook.add_chart({'type': 'line'})
            # Serie 1: Demanda Hogar (Col N)
            chart2.add_series({
                'name':           f"='{sheet_name}'!{col_user_hogar_D}$1", # 'D_Hogar_i{i}(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                'values': f"='{sheet_name}'!{col_user_hogar_D}$2:{col_user_hogar_D}${data_rows}",
                'line':       {'color': 'blue'},
            })
            # Serie 2: Demanda Vehiculo (Col M)
            chart2.add_series({
                'name':           f"='{sheet_name}'!{col_user_sv_D}$1", # 'D_sv_i{i}(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                'values': f"='{sheet_name}'!{col_user_sv_D}$2:{col_user_sv_D}${data_rows}",
                'line':       {'color': 'green'},
            })
            chart2.set_title({'name': f'User {i} Demand (Home vs. Vehicle)'})
            chart2.set_x_axis({'name': 'Period (hour)'})
            chart2.set_y_axis({'name': 'Demand (kWh)'})
            worksheet_horarios.insert_chart('U18', chart2, {'x_scale': chart_x_scale, 'y_scale': chart_y_scale})

            # =================================================
            # GRAFICO 3 - Intercambio Red
            # =================================================
            chart3 = workbook.add_chart({'type': 'line'})
            # Serie 1: Compra Red (Col C)
            chart3.add_series({
                'name':           f"='{sheet_name}'!{col_rh_comunal}$1", # 'rh_Comunal(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                'values': f"='{sheet_name}'!{col_rh_comunal}$2:{col_rh_comunal}${data_rows}",
                'line':       {'color': 'red'},
            })
            # Serie 2: Venta Red (Col D)
            chart3.add_series({
                'name':           f"='{sheet_name}'!{col_hr_comunal}$1", # 'hr_Comunal(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                'values': f"='{sheet_name}'!{col_hr_comunal}$2:{col_hr_comunal}${data_rows}",
                'line':       {'color': 'green'},
            })
            chart3.set_title({'name': 'Community Grid Exchange (Bought vs. Sold)'})
            chart3.set_x_axis({'name': 'Period (hour)'})
            chart3.set_y_axis({'name': 'Energy (kWh)'})
            worksheet_horarios.insert_chart('AD2', chart3, {'x_scale': chart_x_scale, 'y_scale': chart_y_scale})

            # =================================================
            # GRAFICO 4 - SOC BESS (Condicional)
            # =================================================
            if kBESS_opt > 0.1: # Solo mostrar si hay BESS
                chart4 = workbook.add_chart({'type': 'line'})
                chart4.add_series({
                         'name':           f"='{sheet_name}'!{col_socBESS}$1", # 'socBESS(kWh)'
                         'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                         'values': f"='{sheet_name}'!{col_socBESS}$2:{col_socBESS}${data_rows}",
                    'line':       {'color': 'orange'},
                })
                chart4.set_title({'name': 'Community BESS State of Charge'})
                chart4.set_x_axis({'name': 'Period (hour)'})
                chart4.set_y_axis({'name': 'Energy (kWh)'})
                worksheet_horarios.insert_chart('AD18', chart4, {'x_scale': chart_x_scale, 'y_scale': chart_y_scale})

            # =================================================
            # GRAFICO 5 - PV vs Demanda Total
            # =================================================
            chart5 = workbook.add_chart({'type': 'line'})
            # Serie 1: PV Gen (Col E)
            chart5.add_series({
                'name':           f"='{sheet_name}'!{col_ph_comunal}$1", # 'ph_Comunal(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                'values': f"='{sheet_name}'!{col_ph_comunal}$2:{col_ph_comunal}${data_rows}",
                'line':       {'color': '#FFC300'}, # Amarillo
            })
            # Serie 2: Demanda Total Hogar (Col K)
            chart5.add_series({
                'name':           f"='{sheet_name}'!{col_total_hogar_D}$1", # 'Total_Hogar_Demand(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                'values': f"='{sheet_name}'!{col_total_hogar_D}$2:{col_total_hogar_D}${data_rows}",
                'line':       {'color': 'purple'},
            })
            # Serie 3: Demanda Total Vehiculo (Col L)
            chart5.add_series({
                'name':           f"='{sheet_name}'!{col_total_sv_D}$1", # 'Total_SV_Demand(kWh)'
                'categories': f"='{sheet_name}'!{col_periodo}$2:{col_periodo}${data_rows}",
                'values': f"='{sheet_name}'!{col_total_sv_D}$2:{col_total_sv_D}${data_rows}",
                'line':       {'color': 'gray'},
            })
            chart5.set_title({'name': 'Community PV Generation vs. Total Demand'})
            chart5.set_x_axis({'name': 'Period (hour)'})
            chart5.set_y_axis({'name': 'Energy (kWh)'})
            worksheet_horarios.insert_chart('AQ2', chart5, {'x_scale': chart_x_scale, 'y_scale': chart_y_scale})


    print("\n===== REPORTE GUARDADO =====")
    print(f"\nArchivo exportado correctamente en:\n{ruta_salida}")
    print(f"Generado el: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")

else:
    print("\nNo se encontro solucion optima.")
    
    if m.status == GRB.INFEASIBLE:
        print("Estado del modelo: INFACTIBLE. Calculando el IIS para encontrar el conflicto...")
        m.computeIIS()
        iis_path = os.path.join(os.path.dirname(ruta_salida), "conflicto_iis.ilp")
        m.write(iis_path)
        
        print("\n--- INICIO REPORTE IIS ---")
        print(f"Se ha guardado un archivo de conflicto en: {iis_path}")
        
        print("\nLas siguientes RESTRICCIONES forman parte del conflicto:")
        for constr in m.getConstrs():
            if constr.IISConstr:
                print(f"  - {constr.ConstrName}")
        
        print("\nLos siguientes LIMITES DE VARIABLES forman parte del conflicto:")
        for var in m.getVars():
            if var.IISLB > 0:
                print(f"  - Limite inferior de '{var.VarName}'")
            if var.IISUB > 0:
                print(f"  - Limite superior de '{var.VarName}'")
        print("--- FIN REPORTE IIS ---")
        print("\nRevisa las restricciones y limites listados arriba para identificar la contradiccion en tu modelo.")

    else:
        print(f"Estado del modelo: {m.status}")