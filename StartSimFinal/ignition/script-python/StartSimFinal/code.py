Reservoar_h_init = 2.8

while system.tag.read("StartSim").value == True:
        
    #Bibliotek
    import math
    import time
    import random
    import array
    import sys
    
    # Oppretter Pumpestyringsprofil array 
    timer = 24
    styringsProfilArray = array.array('d', [0] * timer)
    styringsProfilSekundforSekund = []
    #Endrer styringsprofilarray fra timer til sekunder
    for i in range(timer):
        styringsProfilArray[i] = system.tag.readBlocking("Pumpestyringsprofil")[0].value[i]
        for t in range(60*60):
            styringsProfilSekundforSekund.append(styringsProfilArray[i])

    #%% Instrument parameter
    MaaleAvvikNivaa = 0.5 # [%]
    MaaleAvvikStromning = 0.5 # [%]
    
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
    #Reservoar_h_init = 2.8 # [m]
    
    Reservoar_lavkutt = 1.4 # [m]
    Reservoar_hoykutt = 5.6 # [m]
    
    #%% Pumpe parameter
    pumpe_q_max = system.tag.read("PumpeMaksStromning").value/ 3600 # 0.033 [m^3/s] - Volumetrisk strømning ved maks pådrag. Funnet ved 120m^3 / 60s*60s
    pumpe_q_min = 0 # [m^3/s] - Volumetrisk strømning ved 4 mA
    pumpe_q_init = 0.015 # [m^3/s] - Volumetrisk strømning ved første iterasjon | funnet ved totalstrømning fra begge pumpene / timer / 60*60
    # 27.01.23 En forenklet versjon som ikke tar for seg mottrykk (net head) i røret og antatt linear pumpekurve
    
    def ReturnK_p(pumpe_q_max, u_range):
        k_p = (pumpe_q_max - pumpe_q_min) / u_range
        return k_p
    
    k_p = ReturnK_p(pumpe_q_max, u_range) # Forholdet mellom pådragssignal og volumetrisk strømning
    
        
    # Filter Parameter
    TidsIntervall = 100
    NivaaVindu = array.array('d', [Reservoar_h_init] * TidsIntervall)
    StromningsVindu = array.array('d', [pumpe_q_init] * TidsIntervall)
        
    
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
    
    #PriReg_settpunkt = 2.8
    PriReg_settpunkt = system.tag.read("Settpunkt").value
    PriReg_settpunkt = PriReg_settpunkt / (100 * Reservoar_h_max)
    
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
    f_maaling_forsinkelse = int(f_forsinkelse/10)
    N_forsinkelse_maaling = int(round(f_maaling_forsinkelse/1)) + 1
    f_maaling_forsinkelse_array = array.array('f', [0] * N_forsinkelse_maaling)
    
    #%% Metoder
    def clip(x, xmin, xmax):
        return max(min(x, xmax), xmin)
    
    def Sum(array):
        sum = 0
        for i in array:
            sum += i
        return sum
        
    # Klasser
    class DataLagring:
        def __init__(self, t, h_t, f_inn_t, f_ut_t, settpunkt):
            self.t = t # Blir brukt som tid og indeks
            self.h_t = h_t # Høyden i bassenget som funksjon av tid
            self.f_inn_t = f_inn_t # Flow in som funksjon av tid, fått av pump.Flow()
            self.f_ut_t = f_ut_t # f_ut_t # flow ut som funksjon av tid, fått av Reservoar.FlowUt()
            self.settpunkt = settpunkt
            self.u = pumpe_q_init 
            
            # Arrays til plotting | simulering
            self.t_array = array.array('f', [0] * t)
            self.h_t_array = array.array('f', [0] * t)
            self.f_inn_t_array = array.array('f', [0] * t)
            self.f_ut_t_array = array.array('f', [0] * t)
            self.settpunkt_array = array.array('f', [0] * t)
            self.u_f_array = array.array('f', [0] * t)
            self.u_h_array = array.array('f', [0] * t)
            self.e_array = array.array('f', [0] * t)
            self.SekReg_u_i_forrige_array = array.array('f', [0] * t)
            self.PriReg_u_i_forrige_array = array.array('f', [0] * t)
        
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
            
    class MiddelVerdiFilter:
        def __init__(self, T_v):
            self.T_v = T_v
        
        def Filtrer(self, vindu):
            summ = Sum(vindu)
            MiddelVerdi = summ / self.T_v
            return MiddelVerdi
               
    #%% Opprett objekter av klassene
    mittReservoar = Reservoar(Reservoar_h_max, Reservoar_h_min, Reservoar_radius, Reservoar_h_init)
    PrimaerRegulator = PrimaerRegulator(PriReg_settpunkt, mittReservoar.h, PriReg_P, PriReg_I, PriReg_D, PriReg_u_man, PriReg_modus)
    SekundaerRegulator = SekundaerRegulator(PrimaerRegulator.Paadrag(), SekReg_pv_init, SekReg_P, SekReg_I, SekReg_D, SekReg_u_man, SekReg_modus)
    minPumpe = Pumpe(k_p, u_min, u_max, SekundaerRegulator.Paadrag())
    minDataLagring = DataLagring(t_stop, mittReservoar.h, minPumpe.FlowInn(), mittReservoar.StatiskForbruk(), PrimaerRegulator.settpunkt)
    minNivaamaaler = Nivaamaaler(mittReservoar.h, MaaleAvvikNivaa)
    minStromningsmaaler = Stromningsmaaler(minPumpe.FlowInn(), MaaleAvvikStromning)
    middelVerdiFilter = MiddelVerdiFilter(TidsIntervall)
    
    #%% Program-loop
    
    # Denne løkken simulerer fullt innløp
    i = 0
    while i < f_maaling_forsinkelse and system.tag.read("StartSim").value == True:
        f_maaling_forsinkelse_array[i] = init
        i += 1
    
    t = t_start
    while t < t_stop and system.tag.read("StartSim").value == True:
            
        # forsinker innstrømning i tanken med 120 iterasjoner | bare valgt et tall her
        PrimaerRegulator.settpunkt = (system.tag.read("Settpunkt").value/100) * 7
        speed = system.tag.read("Speed").value
        f_inn_forsinket = forsinkelse_array[-1]
        forsinkelse_array[1:] = forsinkelse_array[0:-1]
        forsinkelse_array[0] = minPumpe.FlowInn()
        
        # dh_dt = f_in - f_out / A
        dh_dt = (f_inn_forsinket - mittReservoar.StatiskForbruk()) / mittReservoar.Areal() # | *60 for minutter istedet for sekunder
        
        # Mitt Reservoar.h oppdaterer høyden i tanken per iterasjon
        mittReservoar.h += dh_dt
        
        # Måler nivå med Nivåmåler
        minNivaamaaler.prosessverdi = mittReservoar.h
        
        # Filtrerer målt Nivåverdi
        filtrertNivaa = MiddelVerdiFilter.Filtrer(middelVerdiFilter,NivaaVindu)
        NivaaVindu[1:] = NivaaVindu[0:-1]
        NivaaVindu[0] = minNivaamaaler.PVstoy()

        # Oppdater outer loop regulatoren PV | med støy
        PrimaerRegulator.prosessverdi = filtrertNivaa
        
        # Oppdater outer loop error
        PrimaerRegulator.e = PrimaerRegulator.Error()
        
        # Send Primaer u til inner loop SP
        SekundaerRegulator.settpunkt = PrimaerRegulator.Paadrag()
        
        # Forsinker strømningsmåleren med 12 sekunder
        f_maaling_forsinket = f_maaling_forsinkelse_array[-1]
        f_maaling_forsinkelse_array[1:] = f_maaling_forsinkelse_array[0:-1]
        f_maaling_forsinkelse_array[0] = minPumpe.FlowInn()
        
        # Måler strømning med strømningsmåler
        minStromningsmaaler.prosessverdi = f_maaling_forsinket
        
        # Filtrerer målt Strømningsverdi
        filtrertStromning = MiddelVerdiFilter.Filtrer(middelVerdiFilter,StromningsVindu)
        StromningsVindu[1:] = StromningsVindu[0:-1]
        StromningsVindu[0] = minStromningsmaaler.PVstoy()
        
        # Oppdater sekundær regulator PV
        SekundaerRegulator.prosessverdi = filtrertStromning
        
        # Oppdater inner loop error
        SekundaerRegulator.e = SekundaerRegulator.Error()
        
        # Generer nytt pådrag til pumpa
        if styringsProfilSekundforSekund[t] == 1:
                minPumpe.u = SekundaerRegulator.Paadrag()
        elif SekundaerRegulator.modus == False:
                minPumpe.u = SekundaerRegulator.Paadrag()
        else:
            if PrimaerRegulator.prosessverdi < Reservoar_lavkutt:
                for t_m in range(t, t+1800): # t + en halvtime
                    styringsProfilSekundforSekund[t_m] = 1
            elif PrimaerRegulator.prosessverdi > Reservoar_hoykutt:
                minPumpe.u = SekundaerRegulator.Paadrag()
            else:
                minPumpe.u = 0
                
       # Leser settpunkt i realtime, gjør det mulig å oppdatere samtidig som simulatoren kjører
        PrimaerRegulator.settpunkt = system.tag.read("Settpunkt").value
             
        # Lagring av verdier til Datalagringsklassen
        minDataLagring.t_array[t] = t
        minDataLagring.h_t_array[t] = PrimaerRegulator.prosessverdi
        minDataLagring.f_inn_t_array[t] = SekundaerRegulator.prosessverdi
        minDataLagring.f_ut_t_array[t] = mittReservoar.StatiskForbruk()
        minDataLagring.settpunkt_array[t] = PrimaerRegulator.settpunkt
        minDataLagring.u_f_array[t] = SekundaerRegulator.Paadrag()
        minDataLagring.u_h_array[t] = PrimaerRegulator.Paadrag()
        minDataLagring.e_array[t] = PrimaerRegulator.e
        minDataLagring.SekReg_u_i_forrige_array[t] = SekundaerRegulator.u_i_forrige
        minDataLagring.PriReg_u_i_forrige_array[t] = PrimaerRegulator.u_i_forrige
        
         #Gjør om Verdier til Prosent
        Nivaa_prosent = (minDataLagring.h_t_array[t]/7)*100
        PumpePaadrag = minPumpe.u/(u_max-u_min)*100
        if styringsProfilSekundforSekund[t] == 0 and PriReg_modus == True and SekReg_modus == True:
            FlowInn = 0
        else:
            FlowInn = minDataLagring.f_inn_t_array[t]*3600
        FlowUt = minDataLagring.f_ut_t_array[t]*3600
        time.sleep(1/speed)
        
        # Skriver til Tag Nivaa
        system.tag.writeBlocking(["Nivaats"],[Nivaa_prosent])
        system.tag.writeBlocking(["PumpePaadrag"],[PumpePaadrag])
        system.tag.writeBlocking(["FlowInn"],[FlowInn])
        system.tag.writeBlocking(["FlowUt"],[FlowUt])
        t += 1

    

while system.tag.read("StartSim").value == False:
    system.tag.writeBlocking(["Nivaats"],(Reservoar_h_init/7)*100)
    system.tag.writeBlocking(["PumpePaadrag"], 0)
    system.tag.writeBlocking(["FlowInn"], 0)
    system.tag.writeBlocking(["FlowUt"], 0)


 
