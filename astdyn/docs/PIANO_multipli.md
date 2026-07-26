# Asteroidi multipli — stato, bersaglio, piano

## 1. Cosa il software permette oggi

**Funziona già**, validato su Haumea con Hi'iaka e Namaka (report del maggio
2026): tre tracce d'ombra distinte, export SVG e KML.

| componente | stato |
|------------|-------|
| `RelativeMultiBodyPropagator` | integra `[r0, v0, rho1, drho1, ...]`: posizione assoluta del primario piu' vettori relativi dei satelliti. E' l'approccio corretto per i binari stretti — integrare la posizione relativa evita la perdita di precisione della differenza fra due numeri grandi. |
| `MultiBodyPropagator` | corpi mutuamente interagenti in coordinate assolute |
| `ChebyshevEphemerisManager::add_system_body` | carica un corpo da SPK e ne costruisce i polinomi di Chebyshev |
| `OccultationLogic::find_system_occultations` | cerca le occultazioni trattando i corpi **come sistema**, non come oggetti indipendenti |
| proprieta' fisiche | diametro e magnitudine per corpo dalla configurazione (`system_bodies`), con avviso quando si usano i segnaposto |

**Il limite.** `add_system_body` accetta **solo** un `SPKReader`: l'unica via
d'ingresso e' un kernel con posizioni gia' calcolate. Quei kernel esistono per
pochi sistemi ben studiati (Plutone, Haumea); per la gran parte dei binari noti
non esistono affatto.

Entrambi i propagatori accettano `MultiBodyState`, cioe' **stati cartesiani**:
nessuno accetta elementi kepleriani.

## 2. Il bersaglio

Predire le occultazioni dei satelliti partendo dai **parametri orbitali
pubblicati** — l'archivio di Johnston raccoglie semiasse dell'orbita mutua,
eccentricita', inclinazione e periodo per centinaia di binari — e non solo dai
kernel SPK.

Configurazione da far funzionare:

```yaml
objects:
  asteroids: "22"

system_bodies:
  "22":
    name: "Kalliope"
    diameter_km: 150
    H: 6.45
  "22-S1":
    name: "Linus"
    diameter_km: 28
    H: 10.4
    orbit:
      a_km: 1099
      e: 0.005
      i_deg: 103.7
      node_deg: 100.6
      peri_deg: 0.0
      M_deg: 0.0
      epoch_mjd: 61200.0
      period_days: 3.5951
```

Passando **da pochi sistemi a centinaia**.

## 3. Cosa manca

Un solo anello, in tre pezzi.

### 3a. Conversione elementi dell'orbita mutua -> stato relativo

Da (a, e, i, nodo, peri, M) dell'orbita del satellite attorno al primario allo
stato cartesiano relativo. E' la conversione kepleriano->cartesiano che la
libreria gia' esegue, applicata al problema a due corpi satellite-primario.

Il punto delicato e' **mu del sistema**. Si ricava dal periodo con la terza
legge di Keplero:

    mu = 4 pi^2 a^3 / P^2

Preferibile a stimarlo dalle masse: il periodo e' misurato con precisione dalle
curve di luce, mentre le masse dei binari sono spesso incerte al 20-30%.

Attenzione al **piano di riferimento**: i parametri pubblicati sono spesso
riferiti all'equatore J2000 o al piano del cielo, non all'eclittica. La
conversione va esplicitata e documentata, altrimenti l'orbita risulta ruotata.

### 3b. Seconda via d'ingresso nel manager

`ChebyshevEphemerisManager` conserva puntatori a `IChebyshevEphemeris`, e
`SPKChebyshevEphemeris` ne e' una implementazione: basta aggiungerne un'altra
che, invece di leggere da SPK, propaghi con `RelativeMultiBodyPropagator` e
produca gli stessi polinomi.

    void add_orbiting_body(
        const std::string& id,
        const std::string& primary_id,
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& mutual_orbit,
        time::EpochTDB start, time::EpochTDB end, int degree = 12);

### 3c. Lettura della configurazione

Il blocco `system_bodies.<id>.orbit`, accanto a `diameter_km` e `H` che gia'
leggiamo.

## 4. Piano, con verifica a ogni passo

### Passo 1 — conversione e mu dal periodo  [piccolo]
Funzione che da (a, P) ricava mu, e da (a, e, i, nodo, peri, M, mu) lo stato
relativo.

*Verifica:* su un'orbita circolare, |rho| = a e |drho| = 2 pi a / P entro il
numerico. Su Linus: a = 1099 km, P = 3.5951 d -> velocita' relativa attesa
~0.0222 km/s.

### Passo 2 — implementazione di IChebyshevEphemeris per orbite propagate  [medio]
Propaga il satellite sull'intervallo e costruisce i polinomi. Il primario deve
essere gia' nel manager.

*Verifica:* rivalutare i polinomi agli istanti di costruzione e confrontare con
la propagazione diretta; lo scarto deve restare sotto il chilometro.

### Passo 3 — lettura della configurazione e innesto  [piccolo]
`add_orbiting_body` quando il corpo ha un blocco `orbit`, `add_system_body`
quando c'e' un kernel.

*Verifica:* una campagna con la configurazione dell'esempio produce due tracce.

### Passo 4 — validazione contro un caso noto  [il passo che conta]
Confronto con una predizione esterna dello stesso evento, oppure con un evento
gia' osservato di cui si conosca l'esito.

*Verifica:* separazione fra le due tracce e istanti di contatto confrontabili
con il riferimento.

> Senza il passo 4 non sappiamo se il risultato sia giusto: e' la lezione delle
> sessioni precedenti, dove quattro difetti reali erano invisibili ai test
> interni ed emersi solo dal confronto con AstDyS e JPL.

### Passo 5 — kernel SPK dove esistono
Il percorso attuale resta, per i sistemi che li hanno: sono piu' accurati di
un'orbita kepleriana propagata, perche' includono le perturbazioni misurate.

## 5. Limiti da documentare

**Un'orbita kepleriana e' un'approssimazione.** Le orbite mutue dei binari sono
perturbate dalla non sfericita' del primario (J2 del corpo, che per un asteroide
irregolare e' grande), dalle maree e dal Sole. Su tempi lunghi il nodo e il
perielio precedono sensibilmente. La predizione e' attendibile vicino all'epoca
degli elementi e degrada allontanandosene — va detto, non scoperto.

**L'incertezza sulla posizione del satellite** e' spesso maggiore di quella
sull'orbita eliocentrica del sistema. L'ellisse di predizione dovrebbe tenerne
conto; oggi non lo fa.
