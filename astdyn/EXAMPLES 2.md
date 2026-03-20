# AstDyn Core: Esempi d'Uso (EXAMPLES)

Questo documento illustra l'utilizzo dei moduli principali della libreria AstDyn C++, dimostrando come la sua architettura *strongly-typed* permetta di risolvere problemi complessi di astrodinamica in modo elegante, compatto e sicuro.

## 1. Setup Temporale e Scale di Tempo

L'oggetto `Instant` è il cuore della gestione temporale in AstDyn. Nasconde la complessità delle scale di tempo fornendo factory methods sicuri ed evitando operazioni su banali `double` non tipizzati.

```cpp
#include <astdyn/time/epoch.hpp>
#include <iostream>

using namespace astdyn::time;

void esempio_setup_temporale() {
    // Creazione di un istante di tempo a partire da una data MJD (UTC) nota
    // Il compiler ci forza ad esplicitare la scala di base tramite il tipo EpochUTC
    EpochUTC time_utc = EpochUTC::from_mjd(60000.0);

    // I calcoli orbitali richiedono sempre il Terrestrial Time (TT).
    // I tipi Epoch incapsulano le operazioni di conversione rendendole esplicite e sicure.
    EpochTT time_tt = time_utc.to_tt();
    double tempo_tt_val = time_tt.mjd();

    // Visualizziamo la differenza, tipicamente di ~69.184 secondi (fissa post-2017)
    // TAI-UTC = 37s, TT-TAI = 32.184s
    double diff_in_days = time_tt.mjd() - time_utc.mjd();
    double diff_in_seconds = diff_in_days * 86400.0;
    
    std::cout << "Differenza TT - UTC: " << diff_in_seconds << " secondi\n";
}
```

## 2. Trasformazione di Frame Type-Safe

AstDyn codifica i reference frame (ITRF, GCRF, etc.) e le unità di misura direttamente nel *type system* a tempo di compilazione. Non è possibile sommare o sottrarre coordinate affette da basi spaziali differenti.

```cpp
#include <astdyn/types/vectors.hpp>
#include <astdyn/core/frame_tags.hpp>
#include <astdyn/core/units.hpp>
#include <astdyn/coordinates/ReferenceFrame.hpp>
#include <astdyn/propagation/OrbitalElements.hpp>

using namespace astdyn;
using namespace astdyn::core;

void esempio_trasformazione_frame() {
    // Vettore posizione radar in ITRF ("Earth-fixed")
    // Questa espressione è sicura, il tipo porta l'identità del frame!
    types::Vector3<ITRF, Meter> radar_pos_itrf(6378137.0, 0.0, 0.0);
    types::Vector3<ITRF, Meter> radar_vel_itrf(0.0, 0.0, 0.0);

    // Definiamo un'epoca, la rotazione dipende dall'earth orientation (Earth Rotation Era / GMST)
    // Si usa EpochUTC per l'input, poi convertita a TT per i calcoli interni se necessario.
    EpochUTC epoch_utc = EpochUTC::from_mjd(60000.0);
    EpochTT epoch_tt = epoch_utc.to_tt();

    // Estraiamo la matrice di rotazione da ITRF a J2000 (ICRS/GCRF analogo)
    // dal namespace ReferenceFrame interpolando per l'epoca corrente
    Eigen::Matrix3d R_itrf_to_gcrf = coordinates::ReferenceFrame::itrf_to_j2000_simple(epoch_tt);

    // Trasformazione del vettore nello spazio GCRF inerziale
    // N.B: Il compilatore bloccherebbe a compile-time un'operazione semanticamente
    // invalida del tipo: radar_pos_itrf + Vettore3<GCRF, Meter>
    types::Vector3<GCRF, Meter> radar_pos_gcrf(R_itrf_to_gcrf * radar_pos_itrf.to_eigen());
    types::Vector3<GCRF, Meter> radar_vel_gcrf(R_itrf_to_gcrf * radar_vel_itrf.to_eigen());

    // Impacchettiamo le coordinate per convertirle in elementi Kepleriani rispetto alla Terra
    propagation::CartesianElements stato_cartesiano;
    stato_cartesiano.epoch = epoch;
    stato_cartesiano.position = radar_pos_gcrf;
    stato_cartesiano.velocity = radar_vel_gcrf;
    stato_cartesiano.gravitational_parameter = constants::GM_EARTH * 1e9; // m^3/s^2

    // Questa funzione processerà i dati restituendo in maniera pulita la struttura dell'Orbita
    propagation::KeplerianElements stato_kepleriano = propagation::cartesian_to_keplerian(stato_cartesiano);
}
```

## 3. I/O e Propagazione (J2)

Un esempio finale che fa coincidere l'astrazione del lettore di file (I/O) con l'engine di propagazione analitico di J2. Legge uno stato, propaga nello spazio e decodifica elegantemente l'output con la gestione corretta di errori.

```cpp
#include <astdyn/io/OELParser.hpp>
#include <astdyn/propagation/j2_propagator.hpp>
#include <iostream>

using namespace astdyn;

void esempio_io_propagazione() {
    // Utilizziamo l'OELParser per leggere ed estrapolare lo stato dal file target
    // Il parsing usa l'implementazione Monadica "std::expected" (o un Result opzionale)
    // che obbliga il programmatore a gestire fluidamente i possibili IOError di lettura.
    auto risultato_stato = io::OELParser::parse_file("satellite.oel");
    
    if (!risultato_stato.has_value()) {
        std::cerr << "Errore Registrato: " << risultato_stato.error() << "\n";
        return;
    }

    // Estraiamo in totale sicurezza il wrapper dello stato di partenza
    // dopo che la routine ha parsato KEP MJD e validato le coordinate
    propagation::KeplerianElements stato_iniziale = risultato_stato.value();

    // Inizializziamo il propagatore J2 che implementa l'AnalyticalPropagator.
    // Il Propagatore aggiornerà in closed-form l'anomalia vera considerando
    // lo schiacciamento del nucleo terrestre rispetto all'equatore (deriva secolare Nodale).
    propagation::J2Propagator propagatore_j2;

    // Definizione del Target della missione: Epoca Base Iniziale + 12 Ore
    double offset_12h_in_days = 0.5;
    time::EpochTDB nuova_epoca = time::EpochTDB::from_mjd(stato_iniziale.epoch.mjd() + offset_12h_in_days);

    // Propagazione al target prescelto
    auto prop_result = propagatore_j2.propagate(stato_iniziale, nuova_epoca);

    if (prop_result.has_value()) {
        auto stato_futuro = propagation::keplerian_to_cartesian(prop_result.value());

        // Sicurezza a Runtime & Compilazione: Grazie al mapping nativo di "Vector3<GCRF, Meter>",
        // l'IntelliSense e il compilatore garantiscono che l'utente stia operando espressivamente in METRI!
        // Inoltre l'architettura ci ha svincolato da buffer arcaici, non c'è ambiguità pos.x vs pos[0].
        std::cout << "Nuova Posizione X [m]: " << stato_futuro.position.x << "\n";
        std::cout << "Nuova Posizione Y [m]: " << stato_futuro.position.y << "\n";
        std::cout << "Nuova Posizione Z [m]: " << stato_futuro.position.z << "\n";
    }
}
```
