# Optimizacijski pristup i primjena BFGS algoritma u rješavanju problema inverzne kinematike robotskog manipulatora

Ovaj repozitorij sadrži objektno orijentiranu implementaciju rješenja problema inverzne kinematike robotskog manipulatora. Broyden–Fletcher–Goldfarb–Shanno (BFGS) optimizacijski algoritam primijenjen je za rješavanje ovog problema. Ovaj algoritam je implementiran u programskom jeziku C++, u obliku S-funkcije koja je kao blok inverzne kinematike primijenjena u simulaciji koja je rađena u Matlab Simulink-u. Rad bloka inverzne kinematike je u simulaciji demonstriran na KUKA iiwa robotskom manipulatoru koji ima 7 stepeni slobode kretanja.

## Sadržaj
- [Uvod](#uvod)
- [Algoritam](#algoritam)
- [Ulazi i izlazi algoritma](#ulazi-i-izlazi-algoritma)
- [Datoteke](#datoteke)
- [Simulacija](#simulacija)
- [Licenca](#licenca)


## Uvod
Svaki zadatak koji robot obavlja može se svesti na ispravno pozicioniranje i orijentaciju vrha manipulatora, dok izvršenje zadatka obavljaju aktuatori koji pokreću zglobove pri tome mijenjajući varijable zglobova. Direktna kinematika robota izražava ovisnost pozicije i orijentacije  vrha manipulatora od  varijabli zglobova. Inverzna kinematika je postupak određivanja varijabli zglobova robotskog manipulatora na osnovu poznavanja željene pozicije i orijentacije vrha manipulatora. Problem je vrlo nelinearan, što otežava tradicionalna analitička rješenja.  
Za 7-DOF (7 stepeni slobode kretanja) redudantni robotski manipulator rješenje ovog problema nije jedinstveno (u smislu postojanja više konfiguracija zgolobova robota koje daju istu poziciju i orijentaciju vrha manipulatora), što predstavlja dodatni izazov. Ovaj projekt koristi BFGS optimizacijski pristup, kako bi iterativno riješio problem inverzne kinematike. Problem inverzne kinematike je tretiran kao problem nelinearne optimizacije bez ograničenja za varijable, čime je omogućena direktna primjena BFGS metode.


## Algoritam
U inverznoj kinematici, cilj je odrediti zglobne varijable robotskog manipulatora koji rezultiraju željenom pozicijom i orijentacijom vrha manipulatora. Ovaj problem se može postaviti kao optimizacijski problem, gdje je cilj minimizirati grešku između željene pozicije (i orijentacije) i stvarne pozicije (i orijentacije) vrha manipulatora.  Ova formulacija omogućava rješavanje problema inverzne kinematike koristeći optimizacijske tehnike poput BFGS algoritma, koji iterativno minimizira grešku i konvergira ka rješenju koje daje konfiguraciju zglobova za željenu poziciju i orijentaciju. U terminima optimizacije greška pozicije i orijentacije se može izraziti kao funkcija cilja (ili kriterij optimalnosti). BFGS algoritam, kao i algoritam jednodimenzionalnog pretraživana kojeg on koristi, su detaljno opisani u fajlu [Documentation](./Documentation.pdf). 

Funkcija cilja izrazava  mjeru odstupanja trenutne pozicije i orijentacije vrha manipulatora od njegove zeljene pozicije i orijentacije.  
Funkcija cilja u sebe ukljucuje:

- **Grešku pozicije**: izražava rastojanje trenutne pozicije vrha manipulator (koji je funkcija zglobnih varijabli) od zeljene pozicije;
- **Grešku orijentacije**: izražava odstupanje trenutne orijentacije vrha manipulatora od željene.Ovo odstupanje se moze izraziti vektorskim dijelom kvaterniona.

Željena pozicija i željena orijentacija vrha manipulatora se dobiju na osnovu matrice homogene transformacije i ova matrica je  ulaz algoritma. Trenutna poz. i orij. se dobiju iz jednacine direktne  kinematike, koja se dobije na osnovu dh parametara robota (koje cine  i izglobne varijable). U funkciji cilja, zglobne varijable predstavljaju problemske varijable algoritma optimizacije a vektor zglobnih varijabli je izlaz algoritma inverzne kinematike.

## Ulazi i izlazi algoritma

- **Ulazi**:  
  - `Pose`: Željena pozicija i orijentacija vrha manipulatora. Data je kao 4x4 matrica homogene transformacije. Ovu matricu generira blok planiranja trajektorije (koji nije obuhvacen ovim radom); 
  - `Init`: Početna aproksimacija zglobnih varijabli za BFGS algoritam. U svakom koraku  simulacije pocetna vrijednost je jednaka vektoru zglobnih varijabli iz prethodnog koraka.  

- **Izlazi**:  
  - `Config`: Vektor zglobnih varijabli. Predstavlja rezultat racunanja implementiranog algoritma inverzne kinematike. To su varijable zglobova koji omogućavaju postizanje željene pozicije i orijentacije vrha manipulatora.  
  - `Criterion Value`: Vrijednost funkcije kriterija u svakoj tački zadane trajektorije. Izrazava mjeru odstupanja postignute pozicije i orijentacije vrha manipulatora od njegove zeljene pozicije i orijentacije.
    

## Datoteke

Algoritam inverzne kinematike je implementiran u obliku S-funkcije koja je kao blok inverzne kinematike primijenjena u simulaciji koja je rađena u Simulink-u.
- U mapi Libraries u bloku S-funkcije su navedene implemetirane biblioteke:
  - `InverseKinematicsAlgorithm.cpp`: U ovom fajlu se nalazi izvorni kod koji predstavlja implementiciju BFGS optimizacijskog algoritma, funkcije cilja čiji optimum nalazi BFGS algoritam, linijskog pretraživanja  i drugih metoda koje su pozivane u navedenim implementacijama i koje su detaljno opisane u komentaru koda koji se nalazi u navedenom fajlu.    
  - `MatrixLibrary.cpp`: U ovom fajlu su implementirane operacije sa matricama i vektorima kao i funkcije za inicijalizaciju istih.
  
- Header datoteke: `InverseKinematicsAlgorithm.h`, `LineSearch.h`, `MatrixLibrary.h`. U ovim fajlovima su deklarisane klase i metode koje su definisne u istoimenim .cpp fajlovima.
- `•	InverseKinematicsCppOOP_wrapper.cpp`:    U ovom fajlu su zadani DH parametri robota, izvršeno je instanciranje i pozvane su metode implementirane u navedenim bibliotekama.
- `Robot_Params_and_Waypoints_Data.m`: Pored učitavanja parametara modela robota korištenog u simulaciji, izvršavanjem ovog fajla učitavaju se koordinate putnih tačaka, orijentacija vrha manipulatora u putnim tačkama, vremena pristizanja vrha manipulatora u putne tačke, komponente vektora brzine i ubrzanja u putnim tačkama.
- `Robot_Simulation_IK.slx`: U ovom fajlu se implementirana simulacija robota u obluku Simulink blok dijagrama. Pored bloka S-funkcije u kojem je implementiran algoritam inverzne kinematike, u blok dijagramu se nalaze blokovi: Polynomial trajectory – koji na osnovu zadanih putnih tačaka generira trajektoriju; Rotation trajectory- koji za zadanu orijentaciju u  putnim tačkama generira orijentaciju vrha manipulatora duž trajektorije; Visualize Robot-koji omogučava 3D vizualizaciju robota i zadane trajektorije.
- `Documentation.pdf`: Ovaj projekat je rađen kao seminarski rad na Elektrotehničkom fakultetu u Sarajevu. U ovom fajlu je dokumentovan ovaj rad. U njemu su opisani tehnički detalji koji se tiču kinematike robota  i algoritma optimizacije. Opisana je simulacija i rezultati simulacije. U komentaru koda se pozivalo na numerisane dijelove ovog dokumenta. Kod koji je dokumentovan u ovom fajlu je dodatno dorađen u smislu njegove preglednosti a komentar koda je opširniji i detaljniji.
Napomenimo da je ovaj projekat rađen u februaru 2024. godine. U implementaciji je korišten prethodno rađen projekat čiji je izvorni kod dokumentovan u ovon fajlu i koji se nalazi u repozitoriju `Inverse_Kinematics_Cpp_Implementation_and_Simulation_7DOF_Robot`
- `trajExampleUtils.slx` i `visualizeRobot.m`: se koriste u vizualizaciji robota. Ovi fajlovi su preuzeti sa: github.com/ mathworks-robotics/trajectory-planning-robot-manipulators


## Simulacija
Da bi koristili ovaj projekt, potrebno je na računaru imati instaliran MathWorks Matlab softverski paket, verziju R2020a ili noviju verziju. Da bi simulirali i demonstrirali rad implemetiranog algoritma inverzne kinematike potrebno je:
1.	U Matlab-u instalirati C++ kompajler; 
2.	Klonirati repozitorij na računar;
3.	Pokrenuti m-fajl `Robot_Params_and_Waypoints_Data.m`koji se nalazi u ovom repozitoriju; 
4.	Dvoklikom na `Robot_Simulation_IK.slx` fajl otvara se Simulink blok dijagram koji predstavlja  simulaciju robota. Dvoklikom na blok S-funkcije (Inverse Kinematics) i klikom na dugme Build, kompajlira se C++ kod implementiran u s-funkciji.
5.	Klikom na *Run* koje se nalazi u alatnoj traci Simulink-a pokreće se simulacija robotskog manipulatora, pri čemu se otvara prozor sa animacijom robota čiji vrh izvršava zahtjevano kretanje (trajektoriju).



## Licenca
Ovaj projekt je licenciran pod MIT licencom. Detalji licence nalaze se u datoteci [LICENSE](./LICENSE).
