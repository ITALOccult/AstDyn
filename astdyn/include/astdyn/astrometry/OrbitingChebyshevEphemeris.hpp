#ifndef ASTDYN_ASTROMETRY_ORBITING_CHEBYSHEV_EPHEMERIS_HPP
#define ASTDYN_ASTROMETRY_ORBITING_CHEBYSHEV_EPHEMERIS_HPP

#include "astdyn/astrometry/IChebyshevEphemeris.hpp"
#include "astdyn/catalog/CatalogTypes.hpp"
#include "astdyn/core/physics_state.hpp"
#include <vector>
#include <string>

namespace astdyn { struct AstDynConfig; }

namespace astdyn::astrometry {

/**
 * @brief Effemeride di un satellite di asteroide, propagato sulla sua orbita mutua.
 *
 * I sistemi binari pubblicati (archivio di Johnston, letteratura) sono descritti
 * dai parametri dell'orbita del satellite ATTORNO AL PRIMARIO, non da posizioni
 * assolute: per la gran parte dei binari noti non esiste alcun kernel SPK.
 *
 * Questa classe deriva la posizione del satellite sommando alla posizione
 * eliocentrica del primario il vettore relativo, ottenuto risolvendo il problema
 * a due corpi sull'orbita mutua.
 *
 * LIMITE DA CONOSCERE: l'orbita mutua e' trattata come kepleriana. Le orbite dei
 * satelliti di asteroidi sono in realta' perturbate dalla non sfericita' del
 * primario (il J2 di un corpo irregolare e' grande), dalle maree e dal Sole; su
 * tempi lunghi nodo e pericentro precedono sensibilmente. La predizione e'
 * attendibile vicino all'epoca degli elementi e degrada allontanandosene. Dove
 * esiste un kernel SPK e' preferibile, perche' incorpora le perturbazioni misurate.
 */
class OrbitingChebyshevEphemeris : public IChebyshevEphemeris {
public:
    /// Piano cui sono riferiti gli angoli dell'orbita mutua.
    enum class PianoRiferimento {
        Eclittico,    ///< eclittica media J2000 (come il resto della libreria)
        Equatoriale   ///< equatore J2000 — la convenzione piu' diffusa in letteratura
    };

    /**
     * @brief Parametri dell'orbita del satellite attorno al primario.
     *
     * Il parametro gravitazionale del sistema si ricava dal PERIODO tramite la
     * terza legge di Keplero, mu = 4 pi^2 a^3 / P^2, e non dalle masse: il
     * periodo si misura bene dalle curve di luce, mentre le masse dei binari
     * sono spesso incerte del 20-30%.
     */
    struct OrbitaMutua {
        double a_km          = 0.0;   ///< semiasse maggiore [km]
        double e             = 0.0;   ///< eccentricita'
        double i_deg         = 0.0;   ///< inclinazione [gradi]
        double node_deg      = 0.0;   ///< longitudine del nodo ascendente [gradi]
        double peri_deg      = 0.0;   ///< argomento del pericentro [gradi]
        double M_deg         = 0.0;   ///< anomalia media all'epoca [gradi]
        double period_days   = 0.0;   ///< periodo orbitale [giorni]
        time::EpochTDB epoch = time::EpochTDB::from_mjd(0.0);
        PianoRiferimento piano = PianoRiferimento::Equatoriale;
    };

    OrbitingChebyshevEphemeris(
        const physics::KeplerianStateTyped<core::ECLIPJ2000>& primary_elements,
        const OrbitaMutua& orbita,
        time::EpochTDB start_time,
        time::EpochTDB end_time,
        const astdyn::AstDynConfig& config,
        int degree = 12);

    std::tuple<double, double, double> evaluate(time::EpochTDB epoch) const override;

    std::pair<std::tuple<double, double, double>, std::tuple<double, double, double>>
    evaluate_full(time::EpochTDB epoch) const override;

    time::EpochTDB start_epoch() const override { return start_epoch_; }
    time::EpochTDB end_epoch()   const override { return end_epoch_; }
    size_t num_segments()        const override { return segments_.size(); }
    const catalog::ChebyshevSegment& get_segment(time::EpochTDB epoch) const override;

    /// Parametro gravitazionale del sistema [km^3/s^2], ricavato dal periodo.
    static double mu_dal_periodo(double a_km, double period_days);

    /// Vettore posizione del satellite relativo al primario [km], frame eclittico J2000.
    static Eigen::Vector3d vettore_relativo(const OrbitaMutua& orbita, time::EpochTDB t);

private:
    std::vector<catalog::ChebyshevSegment> segments_;
    time::EpochTDB start_epoch_;
    time::EpochTDB end_epoch_;
    double jd_start_ = 0.0;
    double jd_end_   = 0.0;
};

} // namespace astdyn::astrometry

#endif // ASTDYN_ASTROMETRY_ORBITING_CHEBYSHEV_EPHEMERIS_HPP
