while system.tag.read("StartStrom").value == True:	
	#Bibliotek
	import math
	import time
	import random
	import array
	import sys
	
	#%% Pris
	
	pris_kaskade = 0
	pris_smartstyring = 0
	
	# Oppretter Pumpestyringsprofil array 
	timer = 24
	StyringsProfilArray = array.array('d', [0] * timer)
	StyringsProfilSekundforSekund = []
	
	# Oppretter array for strømpriser
	StromPrisArray = array.array('d', [0] * timer)
	StromPrisSekundforSekund = []
	
	#Endrer styringsprofilarray fra timer til sekunder
	for i in range(timer):
	    StyringsProfilArray[i] = system.tag.readBlocking("Pumpestyringsprofil")[0].value[i]
	    StromPrisArray[i] = system.tag.readBlocking("DagensPriser")[0].value[i]
	    for t in range(60*60):
	        StyringsProfilSekundforSekund.append(StyringsProfilArray[i])
	        StromPrisSekundforSekund.append(StromPrisArray[i]) # [Kr/kWt]
	
	#%% Instrument parameter
	MaaleAvvikNivaa = 0.0 # [%]
	MaaleAvvikStromning = 0.0 # [%]
	
	#%% Simulation Time Settings
	t_start = 0  # [s]
	t_stop = 60*60*timer # [s]
	
	#%% Signal parameter
	u_min = 0 # [mA]
	u_max = 16 # [mA]
	u_range = u_max - u_min # [mA]
	
	#%% Reservoar parameter
	Reservoar_h_max = system.tag.read("BassengHoyde").value # 7 [m]
	Reservoar_h_min = 0 # [m]
	Reservoar_radius = system.tag.read("BassengRadius").value # 11 [m]
	Reservoar_h_init = 2.8 # [m]
	
	#%% Pumpe parameter
	pumpe_energiforbruk = system.tag.read("PumpeEffekt").value # kW/t ved max
	
	# kW/h per mA (styresignal)
	pumpe_energi_forbruk_per_mA = pumpe_energiforbruk / 16.0 # [kWt/mA]
	
	# kW/s per mA (Styresignal)
	pumpe_energi_forbruk_per_s_per_mA = pumpe_energi_forbruk_per_mA / 3600 # [kWs/mA]
	
	pumpe_q_max_t = system.tag.read("PumpeMaksStromning").value # 120 [m^3/t] - Volumetrisk strømning ved maks pådrag. Funnet ved 120m^3 / 60s*60s
	pumpe_q_max_s = pumpe_q_max_t / (60*60) # 0.033 [m^3/s]
	pumpe_q_min = 0 # [m^3/s] - Volumetrisk strømning ved 4 mA (0 mA per nå)
	pumpe_q_init = 0.015 # [m^3/s] - Volumetrisk strømning ved første iterasjon | funnet ved totalstrømning fra begge pumpene / timer / 60*60
	# 27.01.23 En forenklet versjon som ikke tar for seg mottrykk (net head) i røret og antatt linear pumpekurve
	
	def ReturnK_p(pumpe_q_max, pumpe_q_min, u_range):
	    k_p = (pumpe_q_max - pumpe_q_min) / u_range
	    return k_p
	
	k_p = ReturnK_p(pumpe_q_max_s, pumpe_q_min ,u_range) # Forholdet mellom pådragssignal og volumetrisk strømning
	
	
	#%% Innløpsventil parameter
	C_v1 = 0.005 # [x] Valve flow coefficient
	
	#%% Utløpsventil parameter
	C_v2 = 0.005 # [x] Valve flow coefficient
	init = round(math.sqrt(Reservoar_h_init) * C_v2, 3)
	
	#%% Primærsløyfe nivå-PI(D)-regulator parameter
	PriReg_P = system.tag.read("ProporsjonalVerdiPrimaer").value
	PriReg_I = system.tag.read("IntegralVerdiPrimaer").value
	PriReg_D = system.tag.read("DerivatVerdiPrimaer").value
	PriReg_u_man = system.tag.read("ManueltPaadragPrimaer").value
	PriReg_settpunkt = 2.8
	PriReg_modus = system.tag.read("PrimaerRegulatorModus").value
	
	
	#%% Sekundærsløyfe strømnings-PI(D)-regulator parameter
	SekReg_P = system.tag.read("ProporsjonalVerdiSekundaer").value
	SekReg_I = system.tag.read("IntegralVerdiSekundaer").value
	SekReg_D = system.tag.read("DerivatVerdiSekundaer").value
	SekReg_u_man = system.tag.read("ManueltPaadragSekundaer").value
	SekReg_modus = system.tag.read("SekundaerRegulatorModus").value
	SekReg_pv_init = pumpe_q_init
	SekReg_sp_init = init
	
	#%% Initialization of time forsinkelse flow in
	f_forsinkelse = 120
	N_forsinkelse_f = int(round(f_forsinkelse/1)) + 1
	forsinkelse_array = array.array('f', [0] * N_forsinkelse_f)
	
	#%% Initialisering av tidsforsinkelse ved måling av strømning
	f_maaling__forsinkelse = int(f_forsinkelse/10)
	N_forsinkelse_maaling_ = int(round(f_maaling__forsinkelse/1)) + 1
	f_maaling__forsinkelse_array = array.array('f', [0] * N_forsinkelse_maaling_)
	
	#%% Klasser 
	def clip(x, xmin, xmax):
	    return max(min(x, xmax), xmin)
	    
	class DataLagring:
	    def __init__(self, t, h_t, f_Inn_t, f_Ut_t, settpunkt):
	        self.t = t # Blir brukt som tid og indeks
	        self.h_t = h_t # Høyden i bassenget som funksjon av tid
	        self.f_Inn_t = f_Inn_t # Flow in som funksjon av tid, fått av pump.Flow()
	        self.f_Ut_t = f_Ut_t # f_Ut_t # flow ut som funksjon av tid, fått av Reservoar.FlowUt()
	        self.settpunkt = settpunkt
	        self.u = pumpe_q_init # Muligens feil
	        
	        # Arrays til plotting | kaskade
	        self.t_array = array.array('f', [0] * t)
	        self.h_t_array = array.array('f', [0] * t)
	        self.f_Inn_t_array = array.array('f', [0] * t)
	        self.f_Ut_t_array = array.array('f', [0] * t)
	        self.settpunkt_array = array.array('f', [0] * t)
	        self.u_f_array = array.array('f', [0] * t)
	        self.u_h_array = array.array('f', [0] * t)
	        self.e_array = array.array('f', [0] * t)
	        self.pris_kaskade_array = array.array('f', [0] * t)
	        self.SekReg_u_i_forrige_array = array.array('f', [0] * t)
	        self.PriReg_u_i_forrige_array = array.array('f', [0] * t)
	        
	        # Arrays til plotting | smartstyring
	        self.h_t_array_smart = array.array('f', [0] * t)
	        self.f_Inn_t_array_smart = array.array('f', [0] * t)
	        self.u_f_array_smart = array.array('f', [0] * t)
	        self.pris_smartstyring_array = array.array('f', [0] * t)
	class Reservoar:
	    def __init__(self, h_max, h_min, radius, h_init):
	        
	        self.h_max = h_max
	        self.h_min = h_min
	        self.radius = radius
	        self.areal = self.Areal()
	        self.h = h_init
	        
	        # Avgrens min og maks nivå i tank
	        self.h = clip(self.h, self.h_min, self.h_max)
	        self.h_prosent = (self.h/self.h_max)*100
	        
	    def Areal(self):
	        
	        Reservoar_areal = math.pi*self.radius**2
	        return Reservoar_areal
	    
	    def StatiskForbruk(self):
	        
	        FlowUt = C_v2 * math.sqrt(self.h)
	        return FlowUt
	        
	class Pumpe:
	    def __init__(self, K_p, u_min, u_max, u):
	        
	        self.K_p = K_p # [%] forsterkning
	        self.u = u # [mA] pådragssignal
	        self.u_max = u_max # [mA]
	        self.u_min = u_min # [mA]
	        
	
	    def FlowInn(self):
	        FlowInn = self.K_p * self.u 
	        return FlowInn # [m^3/s]
	
	class PrimaerRegulator:
	    def __init__(self, settpunkt, prosessverdi, P, I, D, u_man, modus):
	        
	        self.settpunkt = settpunkt
	        self.prosessverdi = prosessverdi
	        self.P = P
	        self.I = I
	        self.D = D
	        self.u_man = u_man
	        self.modus = modus 
	        self.u_i_forrige = u_man 
	        self.e = self.Error()
	        self.u = self.Paadrag()
	        
	    def Error(self):
	        e = self.settpunkt - self.prosessverdi
	        return e
	    
	    def Paadrag(self):
	        if self.modus == True:
	            # Bidrag fra proporsjonalleddet
	            u_p = self.P * self.e
	            
	            # Bidraget fra integralleddet
	            u_i = ((self.P / self.I) * self.e) + self.u_i_forrige
	            
	            # Regn ut total total forsterkning
	            u = u_p + u_i # u_d | ikke implementert
	            
	            #Oppdater u_i_forrige + antiwindup
	            if u > 0.033: # 16
	                self.u_i_forrige = self.u_i_forrige
	            elif u < 0:
	                self.u_i_forrige = self.u_i_forrige
	            else:   
	                self.u_i_forrige = u_i
	            
	        elif self.modus == False:
	            # Manuell drift
	            u = self.u_man
	            
	        else:
	            # Feilmelding
	            pass
	            
	        # Begrens pådragssignal til gitte rammer
	        u = clip(u, 0, 0.033) # u_min, u_max
	        return u
	
	class SekundaerRegulator:
	    def __init__(self, settpunkt, prosessverdi, P, I, D, u_man, modus):
	        
	        self.settpunkt = settpunkt
	        self.prosessverdi = prosessverdi
	        self.P = P
	        self.I = I
	        self.D = D
	        self.u_man = u_man
	        self.modus = modus
	        self.u_i_forrige = u_man 
	        self.e = self.Error()
	        self.u = self.Paadrag()
	        
	    def Error(self):
	        e = self.settpunkt - self.prosessverdi
	        return e
	        
	    def Paadrag(self):
	        if self.modus == True:
	                
	            # Bidrag fra proporsjonalleddet
	            u_p = self.P * self.e
	                
	            # Bidraget fra integralleddet
	            u_i = (self.P / self.I) * self.e + self.u_i_forrige
	                
	            # Regn ut totalforsterkning
	            u = u_p + u_i # + u_d | u_d ikke implementert
	                
	            #Oppdater u_i_forrige + antiwindup
	            if u > 16:
	                self.u_i_forrige = self.u_i_forrige
	            elif u < 0:
	                self.u_i_forrige = self.u_i_forrige
	            else:   
	                self.u_i_forrige = u_i
	            
	        elif self.modus == False:
	            # Manuell drift
	            u = self.u_man
	            
	        else:
	            # Feilmelding
	            pass
	                
	        # Begrens pådrag til gitte rammer
	        u = clip(u, u_min, u_max)
	        return u
	                
	
	class Instrument:
	    def __init__(self, prosessverdi, maaleavvik):
	        
	        self.prosessverdi = prosessverdi
	        self.maaleavvik = maaleavvik
	        
	        def PVstoy(self):
	            pass
	
	class Nivaamaaler(Instrument):
	    def PVstoy(self):
	        PVavvik = (self.prosessverdi * self.maaleavvik) / 100 
	        PVstoy = random.uniform(self.prosessverdi - PVavvik, self.prosessverdi + PVavvik)
	        return PVstoy
	
	class Stromningsmaaler(Instrument):
	    def PVstoy(self):
	        PVavvik = (self.prosessverdi * self.maaleavvik) / 100 
	        PVstoy = random.uniform(self.prosessverdi - PVavvik, self.prosessverdi + PVavvik)
	        return PVstoy
	           
	#%% Opprett objekter av klassene
	mittReservoar = Reservoar(Reservoar_h_max, Reservoar_h_min, Reservoar_radius, Reservoar_h_init)
	PrimaerRegulator = PrimaerRegulator(PriReg_settpunkt, mittReservoar.h, PriReg_P, PriReg_I, PriReg_D, PriReg_u_man, PriReg_modus)
	SekundaerRegulator = SekundaerRegulator(PrimaerRegulator.Paadrag(), SekReg_pv_init, SekReg_P, SekReg_I, SekReg_D, SekReg_u_man, SekReg_modus)
	minPumpe = Pumpe(k_p, u_min, u_max, SekundaerRegulator.Paadrag())
	minDataLagring = DataLagring(t_stop, mittReservoar.h, minPumpe.FlowInn(), mittReservoar.StatiskForbruk(), PrimaerRegulator.settpunkt)
	minNivaamaaler = Nivaamaaler(mittReservoar.h, MaaleAvvikNivaa)
	minStromningsmaaler = Stromningsmaaler(minPumpe.FlowInn(), MaaleAvvikStromning)
	
	
	#%% Program-loop
	
	# Denne løkken simulerer fullt innløp
	for h in range (f_maaling__forsinkelse + 1):
	    f_maaling__forsinkelse_array[h] = init
	
	# Simulerer null forsinkelse på måling i startfasen
	for j in range (f_maaling__forsinkelse + 1):
	    f_maaling__forsinkelse_array[j] = init
	
	
	for sim in range(2):
	    mittReservoar.h = Reservoar_h_init
	    for t in range (t_start, t_stop):
	    
	    	PrimaerRegulator.settpunkt = (system.tag.read("Settpunkt").value/100) * 7
	        # forsinker innstrømning i tanken med 120 iterasjoner | bare valgt et tall her
	        f_Inn_forsinket = forsinkelse_array[-1]
	        forsinkelse_array[1:] = forsinkelse_array[0:-1]
	        forsinkelse_array[0] = minPumpe.FlowInn()
	        
	        # dh_dt = f_in - f_out / A
	        dh_dt = (f_Inn_forsinket - mittReservoar.StatiskForbruk()) / mittReservoar.Areal() # | *60 for minutter istedet for sekunder
	        
	        # Mitt Reservoar.h oppdaterer høyden i tanken per iterasjon
	        mittReservoar.h += dh_dt
	        
	        # Måler nivå med Nivåmåler
	        minNivaamaaler.prosessverdi = mittReservoar.h
	        
	        # Oppdater outer loop regulatoren PV | med støy
	        PrimaerRegulator.prosessverdi = minNivaamaaler.PVstoy()
	        
	        # Oppdater outer loop error
	        PrimaerRegulator.e = PrimaerRegulator.Error()
	        
	        # Send Primaer u til inner loop SP
	        SekundaerRegulator.settpunkt = PrimaerRegulator.Paadrag()
	        
	        # Forsinker strømningsmåling med 12 sekunder
	        f_maaling__forsinket = f_maaling__forsinkelse_array[-1]
	        f_maaling__forsinkelse_array[1:] = f_maaling__forsinkelse_array[0:-1]
	        f_maaling__forsinkelse_array[0] = minPumpe.FlowInn()
	        
	        # Måler strømning med strømningsmåler
	        minStromningsmaaler.prosessverdi = f_maaling__forsinket
	        
	        # Oppdater sekundær regulator PV
	        SekundaerRegulator.prosessverdi = minStromningsmaaler.PVstoy()
	        
	        # Oppdater inner loop error
	        SekundaerRegulator.e = SekundaerRegulator.Error()
	        
	        # Simulerer kaskaderegulering uten smartstyring
	        if sim == 0:
	            
	            # Generer nytt pådrag til pumpa
	            minPumpe.u = SekundaerRegulator.Paadrag()
	            
	            # Lagring av verdier fra kaskaderegulering
	            minDataLagring.h_t_array[t] = minNivaamaaler.PVstoy()
	            minDataLagring.f_Inn_t_array[t] = minStromningsmaaler.PVstoy()
	            minDataLagring.u_f_array[t] = minPumpe.u
	            minDataLagring.pris_kaskade_array[t] = minPumpe.u * pumpe_energi_forbruk_per_s_per_mA * StromPrisSekundforSekund[t]
	            
	            # Lagrer pris
	            pris_kaskade += minDataLagring.pris_kaskade_array[t]
	        
	        # Simulerer kaskaderegulering med smartstyring
	        elif sim == 1:
	            
	            # Generer nytt pådrag til pumpa dersom smartstyringen tillater det
	            if StyringsProfilSekundforSekund[t] == 1:
	                minPumpe.u = SekundaerRegulator.Paadrag()
	             
	            elif SekundaerRegulator == False:
	                minPumpe.u = SekundaerRegulator.Paadrag()
	            else:    
	                minPumpe.u = 0
	                
	            # Lagring av verdier fra smartstyring  
	        minDataLagring.h_t_array_smart[t] = minNivaamaaler.PVstoy()
	        minDataLagring.f_Inn_t_array_smart[t] = minStromningsmaaler.PVstoy()
	        minDataLagring.u_f_array_smart[t] = minPumpe.u
	        minDataLagring.pris_smartstyring_array[t] = minPumpe.u * pumpe_energi_forbruk_per_s_per_mA * StromPrisSekundforSekund[t] # [mA * kWs/mA * Kr/kWs = Kr]
	        #pris_smartstyring += minDataLagring.pris_smartstyring_array[t]   
	        
	
	
	        # Lagring av verdier til Datalagringsklassen
	        minDataLagring.t_array[t] = t
	        #minDataLagring.h_t_array[t] = minNivåmåler.PVstøy()
	        #minDataLagring.f_Inn_t_array[t] = minStrømningsmåler.PVstøy()
	        minDataLagring.f_Ut_t_array[t] = mittReservoar.StatiskForbruk()
	        minDataLagring.settpunkt_array[t] = PrimaerRegulator.settpunkt
	        #minDataLagring.u_f_array[t] = SekundaerRegulator.Pådrag()
	        #minDataLagring.u_f_array[t] = minPumpe.u
	        minDataLagring.u_h_array[t] = PrimaerRegulator.Paadrag()
	        minDataLagring.e_array[t] = PrimaerRegulator.e
	        minDataLagring.SekReg_u_i_forrige_array[t] = SekundaerRegulator.u_i_forrige
	        minDataLagring.PriReg_u_i_forrige_array[t] = PrimaerRegulator.u_i_forrige
	
	# Lagrer pris dette kunne ikke stå i løkken over!!                        
	for t in range (t_start, t_stop):
	    pris_smartstyring += minDataLagring.pris_smartstyring_array[t]
	system.tag.writeBlocking(["Nivaa"],[minDataLagring.h_t_array_smart])
	besparelse = pris_kaskade - pris_smartstyring
	system.tag.writeBlocking(["KaskadePris"], [pris_kaskade])
	system.tag.writeBlocking(["SmartPris"], [pris_smartstyring])
	system.tag.writeBlocking(["Besparelse"], [besparelse])
	system.tag.writeBlocking(["BesparelseProsent"], [(besparelse/pris_kaskade)*100])
	
	    
