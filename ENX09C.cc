// ENX09C.cc — фізична модель детектора згідно експериментальної задачі
//
// Детектор:
//   - Корпус: поліетилен, прямокутний паралелепіпед 53×53×16,3 мм
//   - Чутливий об'єм (сцинтилятор): циліндр Ø39,30×8,55 мм, центр співпадає з центром куба, заповнений порошком
//   - Масова маса порошку: 31,64 мг
//   - Хімічна формула порошку: CoCO3*Co(OH)2*H2O (n=1)
//
// Джерело:
//   - Енергія лінії 60mCo: 58,59 кеВ
//   - Енергія ліній 60Co: 1173,2 кеВ, 1332,5 кеВ
//
// Мета: моделювання проходження гамма-квантів та вторинних частинок через фізичну мішень, збір спектру та просторового розподілу у відповідності до реальної конструкції. Чутливим є лише порошок (сцинтилятор).

#include "G4RunManagerFactory.hh"
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"

#include "G4VUserDetectorConstruction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4VUserActionInitialization.hh"
#include "G4VUserPrimaryGeneratorAction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Gamma.hh"
#include "G4Neutron.hh"
#include "G4GenericIon.hh"
#include "G4IonTable.hh"

#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Event.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4Ions.hh"

#include "G4EmStandardPhysics_option4.hh"
#include "G4EmExtraPhysics.hh"
#include "G4OpticalPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4IonPhysics.hh"


#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <mutex>
#include <filesystem>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <limits>
#include <thread>
#include <cmath>

#include "G4Run.hh"
#include "G4Threading.hh"
#include <atomic>

namespace ENX09C
{
    // Параметри
    inline double BIN_WIDTH = 0.01;   // MeV (дефолт 10 keV)
    inline double MAX_EN = 10.0;   // MeV
    inline int    NUM_BINS = int(std::ceil(MAX_EN / BIN_WIDTH));
    inline double PRIMARY_E = 8.7;    // MeV (джерело можеш задавати /enx/gun/energy)
    inline void Recompute() { NUM_BINS = std::max(1, int(std::ceil(MAX_EN / BIN_WIDTH))); }

    // ADD: форматер десяткової коми (без CSV-розділення)
    static inline std::string fmt_comma(double v, int prec = 6) {
        std::ostringstream oss;
        oss.setf(std::ios::fixed);
        oss << std::setprecision(prec) << v;
        std::string s = oss.str();
        std::replace(s.begin(), s.end(), '.', ',');
        return s;
    }
    // UI-команди для керування бінінгом з макро
    class ParamsMessenger : public G4UImessenger {
    public:
        ParamsMessenger() {
            binWidthCmd = new G4UIcmdWithADoubleAndUnit("/enx/bins/binWidth", this);
            binWidthCmd->SetGuidance("Set histogram bin width");
            binWidthCmd->SetParameterName("dE", false);
            binWidthCmd->SetUnitCategory("Energy");
            binWidthCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

            maxEnergyCmd = new G4UIcmdWithADoubleAndUnit("/enx/bins/maxEnergy", this);
            maxEnergyCmd->SetGuidance("Set histogram max energy");
            maxEnergyCmd->SetParameterName("Emax", false);
            maxEnergyCmd->SetUnitCategory("Energy");
            maxEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

            gammaEnergyCmd = new G4UIcmdWithADoubleAndUnit("/enx/gun/energy", this);
            gammaEnergyCmd->SetGuidance("Set gamma energy");
            gammaEnergyCmd->SetParameterName("Egamma", false);
            gammaEnergyCmd->SetUnitCategory("Energy");
            gammaEnergyCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
        }
        ~ParamsMessenger() override {
            delete binWidthCmd;
            delete maxEnergyCmd;
            delete gammaEnergyCmd;
        }
        void SetNewValue(G4UIcommand* cmd, G4String val) override {
            if (cmd == binWidthCmd) {
                ENX09C::BIN_WIDTH = binWidthCmd->GetNewDoubleValue(val) / CLHEP::MeV;
                ENX09C::Recompute();
            }
            else if (cmd == maxEnergyCmd) {
                ENX09C::MAX_EN = maxEnergyCmd->GetNewDoubleValue(val) / CLHEP::MeV;
                ENX09C::Recompute();
            }
            else if (cmd == gammaEnergyCmd) {
                ENX09C::PRIMARY_E = gammaEnergyCmd->GetNewDoubleValue(val) / CLHEP::MeV;
            }
        }
    private:
        G4UIcmdWithADoubleAndUnit* binWidthCmd{};
        G4UIcmdWithADoubleAndUnit* maxEnergyCmd{};
        G4UIcmdWithADoubleAndUnit* gammaEnergyCmd{};
        G4UIcmdWithABool* writeEntryCSVCmd{};
    };
    inline std::unique_ptr<ParamsMessenger> gParamsMessenger;
    // Анонімний простір імен для глобальних об'єктів, видимих лише у цьому файлі
    namespace {
        // Глобальний потік для запису та м'ютекс для синхронізації доступу
        std::ofstream g_hits_file;
        std::mutex g_hits_file_mutex;
    }

    // --- Чутливий детектор ---
    // SD that counts Co-60 and Co-60m entries into the target (no coordinates)
    class ActivationSD : public G4VSensitiveDetector {
    public:
        ActivationSD(const G4String& name) : G4VSensitiveDetector(name) {}
        G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
            // Count only on entry into the volume to avoid multiple counts per track
            if (step->GetPreStepPoint()->GetStepStatus() != fGeomBoundary) return false;
            auto* def = step->GetTrack()->GetParticleDefinition();
            auto* ion = dynamic_cast<const G4Ions*>(def);
            if (!ion) return false;
            if (ion->GetAtomicNumber() == 27 && ion->GetAtomicMass() == 60) {
                // Any excitation energy > 0 treated as metastable (Co-60m)
                if (ion->GetExcitationEnergy() > 0.0) {
                    co60m_++;
                } else {
                    co60_++;
                }
                return true;
            }
            return false;
        }
        std::pair<int,int> GetCounts() const { return {co60_, co60m_}; }
    private:
        int co60_ = 0;
        int co60m_ = 0;
    };

    // SD that builds a neutron energy histogram at an exit plane
    class NeutronExitSD : public G4VSensitiveDetector {
    public:
        NeutronExitSD(const G4String& name) : G4VSensitiveDetector(name), hist_(NUM_BINS, 0) {}
        G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
            auto* def = step->GetTrack()->GetParticleDefinition();
            if (def != G4Neutron::Neutron()) return false;
            double eMeV = step->GetPreStepPoint()->GetKineticEnergy()/MeV;
            if (eMeV < 0) return false;
            int bin = static_cast<int>(eMeV / BIN_WIDTH);
            if (bin >= 0 && bin < static_cast<int>(hist_.size())) hist_[bin]++;
            return true;
        }
        const std::vector<int>& GetHistogram() const { return hist_; }
    private:
        std::vector<int> hist_;
    };

    // --- Детекторна геометрія ---
    class MyDetector : public G4VUserDetectorConstruction {
    public:
        MyDetector() = default;
        ~MyDetector() override = default;

        G4VPhysicalVolume* Construct() override {
            G4cout << "[MyDetector::Construct] Starting geometry construction." << G4endl;
            auto* nist = G4NistManager::Instance();

            // --- Матеріали ---
            auto* world_mat = nist->FindOrBuildMaterial("G4_Galactic");
            auto* vacuum = world_mat;
            auto* beryllium = nist->FindOrBuildMaterial("G4_Be");
            auto* polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");
            // Додаємо матеріал порошку CoCO3*Co(OH)2*nH2O (приблизна формула)
            // Молярна маса ~ 209.9 г/моль (CoCO3) + 92.94 г/моль (Co(OH)2) + n*18 г/моль (H2O)
            // Для простоти: беремо CoCO3 + Co(OH)2 + 2H2O = 320.84 г/моль
            // Густина: маса/об'єм = 31.64 мг / (π*(19.65 мм)^2*8.55 мм) [мм^3]
            const G4double powder_mass_mg = 31.64; // мг
            const G4double powder_mass_g = powder_mass_mg * 1e-3; // г
            const G4double holeR = 19.65 * mm; // Ø39.3 мм
            const G4double holeH = 8.55 * mm; // повна висота
            constexpr G4double pi_const = 3.14159265358979323846;
            const G4double powder_vol_mm3 = pi_const * holeR * holeR * holeH; // мм^3
            const G4double powder_vol_cm3 = powder_vol_mm3 * 1e-3; // см^3
            const G4double powder_density = powder_mass_g / powder_vol_cm3; // г/см^3
            // Порошок CoCO3*Co(OH)2*H2O
            auto* elCo = nist->FindOrBuildElement("Co");
            auto* elC  = nist->FindOrBuildElement("C");
            auto* elO  = nist->FindOrBuildElement("O");
            auto* elH  = nist->FindOrBuildElement("H");

            auto* powder = new G4Material("CoPowder", 1.5 * g/cm3, 4);
            powder->AddElement(elCo, 2);
            powder->AddElement(elC, 1);
            powder->AddElement(elO, 6);
            powder->AddElement(elH, 4);

            // --- Світ ---
            auto* worldSolid = new G4Box("World", 1 * m, 1 * m, 1 * m);
            auto* worldLV = new G4LogicalVolume(worldSolid, world_mat, "WorldLV");
            auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV", nullptr, false, 0, true);

            // Допоміжна функція для атрибутів візуалізації
            auto wireOnly = [](const G4Colour& c) {
                auto* vis = new G4VisAttributes(c);
                vis->SetForceWireframe(true);
                vis->SetForceSolid(false);
                vis->SetForceAuxEdgeVisible(true);
                return vis;
            };

            // Маркер джерела (залишаємо для візуалізації)
            {
                auto* gunSolid = new G4Sphere("GunSolid", 0, 1.5 * mm, 0, 360 * deg, 0, 180 * deg);
                auto* gunLV = new G4LogicalVolume(gunSolid, vacuum, "GunLV");
                gunLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.7, 1.0)));
                new G4PVPlacement(nullptr, { 0,0,-5 * cm }, gunLV, "GunPV", worldLV, false, 0, true);
            }
            // Повертаємо два послідовних Be-циліндра перед детекторами
            const G4double beR = 5 * cm;
            const G4double beThick = 7 * cm;
            const G4double beZ1 = 2 * cm + beThick / 2.0;
            auto* beSolid1 = new G4Tubs("BeCore1", 0, beR, beThick / 2.0, 0, 360 * deg);
            auto* beLV1 = new G4LogicalVolume(beSolid1, beryllium, "BeLV1");
            beLV1->SetVisAttributes(wireOnly(G4Colour(0.2, 0.6, 0.8)));
            new G4PVPlacement(nullptr, { 0,0,beZ1 }, beLV1, "BePV1", worldLV, false, 0, true);
            // Другий Be-циліндр (щільно за першим, без проміжку)
            const G4double gap = 0.0 * cm; // мінімальний проміжок
            const G4double beZ2 = beZ1 + beThick;
            auto* beSolid2 = new G4Tubs("BeCore2", 0, beR, beThick / 2.0, 0, 360 * deg);
            auto* beLV2 = new G4LogicalVolume(beSolid2, beryllium, "BeLV2");
            beLV2->SetVisAttributes(wireOnly(G4Colour(0.1, 0.4, 0.7)));
            new G4PVPlacement(nullptr, { 0,0,beZ2 }, beLV2, "BePV2", worldLV, false, 0, true);

            // --- Детектор (корпус) ---
            // Створюємо суцільний блок для корпусу
            const G4double detX = 53.0 * mm;
            const G4double detY = 53.0 * mm;
            const G4double detZ = 16.3 * mm;
            auto* detSolidBox = new G4Box("DetectorBox", detX / 2, detY / 2, detZ / 2);
            auto* detLV = new G4LogicalVolume(detSolidBox, polyethylene, "DetectorLV");
            detLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.8, 0.8, 0.8, 0.3)));
            new G4PVPlacement(nullptr, {0, 0, 0}, detLV, "DetectorPV", worldLV, false, 0, true);

            // --- Чутливий об'єм (порошок, розміщений всередині корпусу) ---
            // Створюємо циліндр для порошку і розміщуємо його як дочірній об'єм всередині корпусу.
            // Geant4 автоматично замінить матеріал корпусу (поліетилен) на матеріал порошку в цій області.
            auto* powderSolid = new G4Tubs("Powder", 0, holeR, holeH / 2, 0, 360 * deg);
            fScintillatorLV = new G4LogicalVolume(powderSolid, powder, "ScintillatorLV");
            fScintillatorLV->SetVisAttributes(new G4VisAttributes(G4Colour(0.3, 0.3, 1.0, 0.7)));
            new G4PVPlacement(nullptr, {0, 0, 0}, fScintillatorLV, "ScintillatorPV", detLV, false, 0, true);

            // --- Площина для спектра нейтронів після Be ---
            const G4double SD_HALF_THICK = 0.25 * mm;
            const G4double sdZ = beZ2 + beThick / 2.0 + SD_HALF_THICK;
            auto* sdSolid = new G4Tubs("NeutronSDSolid", 0, beR, SD_HALF_THICK, 0, 360 * deg);
            fNeutronExitLV = new G4LogicalVolume(sdSolid, vacuum, "NeutronExitLV");
            fNeutronExitLV->SetVisAttributes(wireOnly(G4Colour(1.0, 0.0, 0.0)));
            new G4PVPlacement(nullptr, { 0,0,sdZ }, fNeutronExitLV, "NeutronExitPV", worldLV, false, 0, true);

            G4cout << "[MyDetector::Construct] Geometry construction finished." << G4endl;
            return worldPV;
        }

        void ConstructSDandField() override {
            G4cout << "[MyDetector::ConstructSDandField] Constructing SD." << G4endl;
            auto* sdm = G4SDManager::GetSDMpointer();
            auto* act = new ActivationSD("ActivationSD");
            sdm->AddNewDetector(act);
            SetSensitiveDetector(fScintillatorLV, act);

            auto* nsd = new NeutronExitSD("NeutronExitSD");
            sdm->AddNewDetector(nsd);
            SetSensitiveDetector(fNeutronExitLV, nsd);
            G4cout << "[MyDetector::ConstructSDandField] SD constructed and assigned." << G4endl;
        }

    private:
        G4LogicalVolume* fScintillatorLV = nullptr;
        G4LogicalVolume* fNeutronExitLV = nullptr;
    };

    // --- Фізичний список ---
    class MyPhysicsList : public G4VModularPhysicsList {
    public:
        MyPhysicsList() {
            SetVerboseLevel(1);
            RegisterPhysics(new G4DecayPhysics());
            RegisterPhysics(new G4RadioactiveDecayPhysics());
            RegisterPhysics(new G4EmStandardPhysics_option4());
            // Enable gamma-nuclear interactions (photo-nuclear, mu-nuclear, etc.)
            RegisterPhysics(new G4EmExtraPhysics());
            RegisterPhysics(new G4OpticalPhysics());
            RegisterPhysics(new G4HadronElasticPhysics());
            RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
            RegisterPhysics(new G4IonPhysics());
        }
        ~MyPhysicsList() override = default;
    };

    // --- Джерело частинок ---
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    public:
                PrimaryGeneratorAction() {
                        gun_ = new G4ParticleGun(1);
                        gun_->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
                        gun_->SetParticleMomentumDirection({0,0,1});
                        gun_->SetParticlePosition({0,0,-5*cm});
                }
                ~PrimaryGeneratorAction() override { delete gun_; }

        // Повертаємо універсальний спектр (як було)
        static double GammaSpectrumPDF(double x) {
            constexpr double xm = 8.7;
            constexpr double s = -1.1647;
            constexpr double b = -0.00487;
            constexpr double x0 = 1.0;
            constexpr double A = 457.89319;
            double num = std::pow(xm - x, s);
            double den = 1.0 + std::pow(b * (x - x0), 2);
            return (x < x0 || x > xm) ? 0.0 : A * num / den;
        }
        static double SampleGammaEnergy() {
            constexpr double xmin = 1.0;
            constexpr double xmax = 8.7;
            constexpr double pdf_max = 500.0;
            while (true) {
                double x = xmin + (xmax - xmin) * G4UniformRand();
                double y = pdf_max * G4UniformRand();
                if (y < GammaSpectrumPDF(x)) return x;
            }
        }

        void GeneratePrimaries(G4Event* evt) override {
            double energy = SampleGammaEnergy();
            gun_->SetParticleEnergy(energy * MeV);
            gun_->GeneratePrimaryVertex(evt);
        }
        void SetEnergy(double en_MeV) { gun_->SetParticleEnergy(en_MeV*MeV); }
    private:
        G4ParticleGun* gun_;
  };
  // Прогрес % — один друк на крок відсотка (глобально-атомарний)
  class ProgressEventAction : public G4UserEventAction {
  public:
      void EndOfEventAction(const G4Event*) override {
          static std::atomic<int> eventCounter{ 0 };
          static std::atomic<int> lastPercent{ -1 };
          int count = ++eventCounter;
          int total = G4RunManager::GetRunManager()->GetNumberOfEventsToBeProcessed();
          if (total <= 0) return;
          int percent = count * 100 / total;

          int expected = lastPercent.load(std::memory_order_relaxed);
          // Піднімаємо lastPercent до percent, друкуємо лише при реальному зростанні
          while (percent > expected) {
              if (lastPercent.compare_exchange_weak(expected, percent,
                  std::memory_order_relaxed,
                  std::memory_order_relaxed)) {
                  bool useStderr = true;
                  if (useStderr)
                      std::cerr << "[Progress] " << percent << "% (" << count << "/" << total << ")\n";
                  else
                      G4cout << "[Progress] " << percent << "% (" << count << "/" << total << ")\n";
                  break;
              }
              // expected оновлено поточним значенням lastPercent — перевіримо знову
          }
      }
  };
    // AnalysisManager: merges worker data and writes two CSVs
    class AnalysisManager : public G4UserRunAction {
    public:
        AnalysisManager() { isMaster_ = !G4Threading::IsWorkerThread(); }

        void EndOfRunAction(const G4Run*) override {
            if (!isMaster_) {
                auto* sdm = G4SDManager::GetSDMpointer();
                if (auto* nsd = dynamic_cast<NeutronExitSD*>(sdm->FindSensitiveDetector("NeutronExitSD"))) {
                    MergeSpectrum(nsd->GetHistogram());
                }
                if (auto* act = dynamic_cast<ActivationSD*>(sdm->FindSensitiveDetector("ActivationSD"))) {
                    MergeActivation(act->GetCounts());
                }
                return;
            }
            // Master: write files
            std::filesystem::create_directories("results/spectrum");
            {
                std::ofstream csv("results/spectrum/neutron_spectrum.csv");
                csv << "E_MeV;count\n";
                for (size_t i = 0; i < spectrum_.size(); ++i) {
                    double e = (i + 0.5) * BIN_WIDTH;
                    csv << fmt_comma(e) << ";" << spectrum_[i] << "\n";
                }
            }
            std::filesystem::create_directories("results/activation");
            {
                std::ofstream csv("results/activation/co60_stats.csv");
                csv << "Co60;Co60m\n" << co60_ << ";" << co60m_ << "\n";
            }
            G4cout << "[Run] Wrote spectrum and activation stats.\n";
        }

        void MergeSpectrum(const std::vector<int>& workerHist) {
            std::lock_guard<std::mutex> lock(mx_);
            if (spectrum_.size() != workerHist.size()) spectrum_.assign(workerHist.size(), 0);
            for (size_t i = 0; i < workerHist.size(); ++i) spectrum_[i] += workerHist[i];
        }
        void MergeActivation(const std::pair<int,int>& counts) {
            std::lock_guard<std::mutex> lock(mx_);
            co60_ += counts.first;
            co60m_ += counts.second;
        }
    private:
        bool isMaster_ = true;
        std::vector<int> spectrum_;
        int co60_ = 0;
        int co60m_ = 0;
        std::mutex mx_;
    };

  // ✅ Тільки AnalysisManager → ActionInitialization
  class ActionInitialization : public G4VUserActionInitialization {
  public:
      void Build() const override {
          if (!ENX09C::gParamsMessenger) ENX09C::gParamsMessenger = std::make_unique<ENX09C::ParamsMessenger>();
          SetUserAction(new PrimaryGeneratorAction());
          SetUserAction(new ProgressEventAction());
          SetUserAction(new AnalysisManager());
      }
      void BuildForMaster() const override {
          if (!ENX09C::gParamsMessenger) ENX09C::gParamsMessenger = std::make_unique<ENX09C::ParamsMessenger>();
          SetUserAction(new AnalysisManager());
      }
  };
} // namespace ENX09C
// ===== ПІДКЛЮЧЕННЯ ДО ЄДИНОГО МОСТУ ENTERPRISE =====
namespace Enterprise { int StartFromBridgeQt(int argc, char** argv); }
int main(int argc, char** argv) {
    G4cout << "===================================================" << G4endl;
    G4cout << " ENX09C Simulation Program Version 3.0" << G4endl;
    G4cout << " == Enterprise project Version 9.3 ==" << G4endl;
    G4cout << " Modeling gamma → (n,...) on two sequential Be blocks + physical detector ENX09C" << G4endl;
    G4cout << " Gamma source, variable energy, neutron spectrum after Be and in detector" << G4endl;
    G4cout << "===================================================" << G4endl;
    // Парсимо вибір MT // Пріоритет: аргументи -> ENV -> дефолт
    std::cout << "[Enterprise] Enter threads (0 = serial+Qt, 1 = MT 1 thread, >1 = MT N threads): ";
    int threads = 0;
    std::cin >> threads;
    std::cin.ignore(std::numeric_limits<std::streamsize>::max(), '\n');

    bool useMT = (threads > 0);
    auto* runManager = G4RunManagerFactory::CreateRunManager(
        useMT ? G4RunManagerType::MT : G4RunManagerType::Serial);
        if (useMT) {
    auto* mt = dynamic_cast<G4MTRunManager*>(runManager);
        if (mt) {
    const int hw = std::max(1u, std::thread::hardware_concurrency());
    int nThreads = std::max(1, threads);
    mt->SetNumberOfThreads(std::min(nThreads, hw));
    G4cout << "[MT] Requested " << nThreads << " threads; using "
                << mt->GetNumberOfThreads() << " (P = " << hw << ")\n";
        } else {
        G4cout << "[MT] dynamic_cast<G4MTRunManager*> failed! Not running in MT mode?\n";
        }
    }
    runManager->SetUserInitialization(new ENX09C::MyDetector());
    runManager->SetUserInitialization(new ENX09C::MyPhysicsList());
    runManager->SetUserInitialization(new ENX09C::ActionInitialization());
    // Ініціалізуєш сцену ДО запуску будь-якої макро або UI
    runManager->Initialize();
    if (threads == 0) {
        Enterprise::StartFromBridgeQt(argc, argv);
    } else {
        auto* UIm = G4UImanager::GetUIpointer();
        UIm->ApplyCommand("/control/execute run_phase.mac");
    }
    std::cout << "[INFO] Program finished. Press Enter to exit...\n";
    std::cin.get();
    delete runManager;
    return 0;
}