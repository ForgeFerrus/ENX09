#include <mutex>
// ENX09.cc — один Be-циліндр + вхідний диск SD, CSV-вивід (;-розділювач, десяткова кома),
// MT-злиття, прогрес у %, макро-керування, односотовий файл.
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VUserPhysicsList.hh"
#include "G4VUserPrimaryGeneratorAction.hh"
#include "G4VUserActionInitialization.hh"
#include "G4UserRunAction.hh"
#include "G4UserEventAction.hh"
#include "G4SDManager.hh"
#include "G4VSensitiveDetector.hh"
#include "G4UImanager.hh"
#include "G4VisExecutive.hh"
#include "G4UIExecutive.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4Sphere.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4SystemOfUnits.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4PhysicsListHelper.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4Threading.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4GenericMessenger.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ParticleChange.hh"
#include "G4Step.hh"
#include "Randomize.hh"

#include <cstdlib>
#include <optional>
#include <fstream>
#include <iostream>
#include <iomanip>
#include <vector>
#include <string>
#include <cmath>
#include <algorithm>
#include <filesystem>
#include <atomic>
#include <memory>
#include <thread>
#include <map>
#include <limits>
#include <sstream>
#include <corecrt_math_defines.h>


namespace ENX09A 
{
   // Параметри
   inline double BIN_WIDTH = 0.01;   // MeV (дефолт 10 keV)
   inline double MAX_EN = 10.0;   // MeV
   inline int    NUM_BINS = int(std::ceil(MAX_EN / BIN_WIDTH));
   inline double PRIMARY_E = 8.7;    // MeV (джерело можеш задавати /gun/energy)
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

            // writeEntryCSVCmd більше не потрібен, CSV завжди пишеться
        }
        ~ParamsMessenger() override {
            delete binWidthCmd;
            delete maxEnergyCmd;
            // writeEntryCSVCmd більше не потрібен
        }
        void SetNewValue(G4UIcommand* cmd, G4String val) override {
            if (cmd == binWidthCmd) {
                ENX09A::BIN_WIDTH = binWidthCmd->GetNewDoubleValue(val) / CLHEP::MeV;
                ENX09A::Recompute();
            }
            else if (cmd == maxEnergyCmd) {
                ENX09A::MAX_EN = maxEnergyCmd->GetNewDoubleValue(val) / CLHEP::MeV;
                ENX09A::Recompute();
            }
            // writeEntryCSVCmd більше не потрібен
        }
    private:
        G4UIcmdWithADoubleAndUnit* binWidthCmd{};
        G4UIcmdWithADoubleAndUnit* maxEnergyCmd{};
        G4UIcmdWithABool* writeEntryCSVCmd{};
    };

    inline std::unique_ptr<ParamsMessenger> gParamsMessenger;
    // + додай прапорець керування (зробимо його статичним)
    // static std::atomic<bool> gWriteEntryCSV{ false }; // прапорець більше не потрібен
    // SD: ловить gamma перед Be; веде гістограму і лог хітів у CSV per-thread
    class GammaEntrySD : public G4VSensitiveDetector {
    public:
        // Додаємо метод для доступу до гістограми (має бути public)
        const std::vector<int>& GetHistogram() const { return hist_; }
        GammaEntrySD(const G4String& name)
            : G4VSensitiveDetector(name), hist_(ENX09A::NUM_BINS, 0) {
            // Завжди створюємо results і відкриваємо CSV
            std::filesystem::create_directories("results");
            static std::mutex mtx;
            std::lock_guard<std::mutex> lock(mtx);
            static bool opened = false;
            if (!opened) {
                out_.open("results/gamma_entry.csv", std::ios::out | std::ios::trunc);
                out_ << "E_MeV;x_mm;y_mm;z_mm;ux;uy;uz;theta_deg;r_mm\n";
                opened = true;
            } else {
                out_.open("results/gamma_entry.csv", std::ios::out | std::ios::app);
            }
        }

        ~GammaEntrySD() override {
            if (out_.is_open()) out_.close();
        }

        void FlushHistogram() {
            // Більше не пишемо спектр у CSV, тільки детальні дані у gamma_entry.csv
        }

        G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
            const auto* trk = step->GetTrack();
            if (trk->GetParticleDefinition() != G4Gamma::Gamma()) return false;

            const auto* pre = step->GetPreStepPoint();
            // Якщо хочеш рахувати тільки перетини диска, розкоментуй:
            // if (pre->GetStepStatus() != fGeomBoundary) return false;

            const double E = pre->GetKineticEnergy() / MeV;
            const int nb = static_cast<int>(hist_.size());
            int bin = static_cast<int>(E / ENX09A::BIN_WIDTH);
            if (bin >= 0 && bin < nb) hist_[bin]++;

            if (out_.is_open()) {
                const auto pos = pre->GetPosition();
                const double x = pos.x() / mm, y = pos.y() / mm, z = pos.z() / mm;
                const auto dir = pre->GetMomentumDirection();
                const double ux = dir.x(), uy = dir.y(), uz = dir.z();
                static constexpr double pi = 3.14159265358979323846;
                const double theta = std::acos(std::clamp(uz, -1.0, 1.0)) * 180.0 / pi;
                const double r = std::hypot(x, y);
                out_ << ENX09A::fmt_comma(E) << ";"
                    << ENX09A::fmt_comma(x) << ";"
                    << ENX09A::fmt_comma(y) << ";"
                    << ENX09A::fmt_comma(z) << ";"
                    << ENX09A::fmt_comma(ux) << ";"
                    << ENX09A::fmt_comma(uy) << ";"
                    << ENX09A::fmt_comma(uz) << ";"
                    << ENX09A::fmt_comma(theta) << ";"
                    << ENX09A::fmt_comma(r) << "\n";
            }
            return true;
        }

    private:
        std::ofstream out_;
        std::vector<int> hist_;
    };

    class MyDetector : public G4VUserDetectorConstruction {
    public:
        G4VPhysicalVolume* Construct() override {
            auto* nist = G4NistManager::Instance();
            auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
            auto* tungsten = nist->FindOrBuildMaterial("G4_W");
            auto* beryllium = nist->FindOrBuildMaterial("G4_Be");

            auto* worldSolid = new G4Box("World", 1 * m, 1 * m, 1 * m);
            auto* worldLV = new G4LogicalVolume(worldSolid, vacuum, "WorldLV");
            auto* worldPV = new G4PVPlacement(nullptr, {}, worldLV, "WorldPV", nullptr, false, 0, true);

            auto wireOnly = [](G4Colour c) {
                auto* vis = new G4VisAttributes(c);
                vis->SetForceWireframe(true);
                vis->SetForceSolid(false);
                vis->SetForceAuxEdgeVisible(true);
                return vis;
                };

            // Маркер джерела
            {
                auto* gunSolid = new G4Sphere("GunSolid", 0, 1.5 * mm, 0, 360 * deg, 0, 180 * deg);
                auto* gunLV = new G4LogicalVolume(gunSolid, vacuum, "GunLV");
                gunLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.7, 1.0)));
                new G4PVPlacement(nullptr, { 0,0,-10 * cm }, gunLV, "GunPV", worldLV, false, 0, true);
            }

            // Вольфрамова мішень (93×55×2 мм) у z = −5 см
            {
                auto* tgtSolid = new G4Box("Target", 93 * mm / 2, 55 * mm / 2, 1 * mm);
                auto* tgtLV = new G4LogicalVolume(tgtSolid, tungsten, "TargetLV");
                new G4PVPlacement(nullptr, { 0,0,-5 * cm }, tgtLV, "TargetPV", worldLV, false, 0, true);
                tgtLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.4, 0.4)));
            }

            // Be-циліндр
            const G4double beR = 5 * cm;
            const G4double beThick = 7 * cm;
            const G4double beZ = 2 * cm + beThick / 2.0;
            auto* beSolid = new G4Tubs("BeCore", 0, beR, beThick / 2.0, 0, 360 * deg);
            auto* beLV = new G4LogicalVolume(beSolid, beryllium, "BeLV");
            beLV->SetVisAttributes(wireOnly(G4Colour(0.2, 0.6, 0.8)));
            new G4PVPlacement(nullptr, { 0,0,beZ }, beLV, "BePV", worldLV, false, 0, true);

            // Вхідний SD-диск (вакуум) перед торцем
            const G4double SD_HALF_THICK = 0.25 * mm;
            const G4double sdZ = beZ - beThick / 2.0 - SD_HALF_THICK;
            auto* sdSolid = new G4Tubs("GammaEntryDisc", 0, beR, SD_HALF_THICK, 0, 360 * deg);
            fGammaEntryLV = new G4LogicalVolume(sdSolid, vacuum, "GammaEntryLV");
            fGammaEntryLV->SetVisAttributes(wireOnly(G4Colour(1.0, 1.0, 0.0)));
            new G4PVPlacement(nullptr, { 0,0,sdZ }, fGammaEntryLV, "GammaEntryPV", worldLV, false, 0, true);

            return worldPV;
        }

        void ConstructSDandField() override {
            auto* SDM = G4SDManager::GetSDMpointer();
            auto* sd = new GammaEntrySD("GammaEntrySD");
            SDM->AddNewDetector(sd);
            fGammaEntryLV->SetSensitiveDetector(sd);
        }

    private:
        G4LogicalVolume* fGammaEntryLV = nullptr;
  };
  // ── MyPhotoNuclearProcess
  // Цей клас реалізує процес фотоядерної взаємодії
  class MyPhotoNuclearProcess : public G4VDiscreteProcess {
  public:

      G4ParticleChange fParticleChange; // для вторинних частинок
      MyPhotoNuclearProcess() : G4VDiscreteProcess("MyPhotoNuclear") {
          SetProcessType(fHadronic);
          SetProcessSubType(131); // fPhotoNuclear
      }
      G4bool IsApplicable(const G4ParticleDefinition* p) { return p == G4Gamma::Gamma(); }
      G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* cond) override {
          *cond = NotForced;
          return 5 * cm; // груба оцінка, або можна зробити енергозалежною
      }

      G4VParticleChange* PostStepDoIt(const G4Track& track, const G4Step&) override {
          fParticleChange.Initialize(track);
          if (track.GetKineticEnergy() < 7.0 * MeV) return &fParticleChange; // ← без крапки з комою

          // Якщо хочеш зберегти гамма-only аналіз, нічого не додаємо.
          // Якщо потрібен тест вторинних — тут додаси нейтрон і повернеш fParticleChange.

          return &fParticleChange;
      }

  private:
      G4ThreeVector RandomDirection() {
          double theta = acos(2 * G4UniformRand() - 1);
          double phi = 2 * M_PI * G4UniformRand();
          return G4ThreeVector(sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta));
      }
  };
  // Фізичний список
  class MyPhysicsList : public G4VUserPhysicsList {
  public:
    void ConstructParticle() override {
      G4Electron::ElectronDefinition();
      G4Gamma::GammaDefinition();
      G4Positron::PositronDefinition();
    }
    void ConstructProcess() override {
      AddTransportation();
      auto* helper = G4PhysicsListHelper::GetPhysicsListHelper();
      auto* eMinus = G4Electron::Electron();
      helper->RegisterProcess(new G4eMultipleScattering(), eMinus);
      helper->RegisterProcess(new G4eIonisation(), eMinus);
      helper->RegisterProcess(new G4eBremsstrahlung(), eMinus);
      auto* gamma = G4Gamma::Gamma();
      helper->RegisterProcess(new G4PhotoElectricEffect(), gamma);
      helper->RegisterProcess(new G4ComptonScattering(),  gamma);
      helper->RegisterProcess(new G4GammaConversion(),    gamma);
      // наш кастомний γ → n процес
      auto* myPhoto = new MyPhotoNuclearProcess(); // з класу вище
      helper->RegisterProcess(myPhoto, gamma);
    }
    void SetCuts() override {
      SetCutValue(0.1*mm, "gamma");
      SetCutValue(0.1*mm, "e-");
      SetCutValue(0.1*mm, "e+");
    }
  };

  class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
  public:
    PrimaryGeneratorAction() {
      gun_ = new G4ParticleGun(1);
      gun_->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("e-"));
      gun_->SetParticleEnergy(PRIMARY_E*MeV);
      gun_->SetParticleMomentumDirection({0,0,1});
      gun_->SetParticlePosition({0,0,-10*cm});
    }
    ~PrimaryGeneratorAction() override { delete gun_; }
    void GeneratePrimaries(G4Event* evt) override { gun_->GeneratePrimaryVertex(evt); }
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
  // ✅ AnalysisManager — єдиний RunAction
  class AnalysisManager : public G4UserRunAction {
  public:
      AnalysisManager() {
          isMaster_ = !G4Threading::IsWorkerThread();
          messenger_ = new G4GenericMessenger(this, "/analysis/", "Analysis control");
          if (isMaster_) {
              messenger_->DeclareMethod("run", &AnalysisManager::RunCSV)
                  .SetGuidance("Merge spectrum + generate CSV").SetParameterName("name", false);
              messenger_->DeclareMethod("exportSpectrumTXT", &AnalysisManager::ExportSpectrumTXT)
                  .SetGuidance("Write text energy spectrum").SetParameterName("name", false);
          }
      }
      ~AnalysisManager() override { delete messenger_; }

      void EndOfRunAction(const G4Run*) override {
          if (!isMaster_) {
              auto* base = G4SDManager::GetSDMpointer()->FindSensitiveDetector("GammaEntrySD");
              if (auto* sd = dynamic_cast<ENX09A::GammaEntrySD*>(base)) {
                  const auto& localHist = sd->GetHistogram();
                  auto* master = const_cast<AnalysisManager*>(static_cast<const AnalysisManager*>(G4MTRunManager::GetMasterRunManager()->GetUserRunAction()));
                  if (master) master->MergeHistogram(localHist);
              }
              return;
          }
          // === Додаємо збереження спектра у файл ===
          std::filesystem::create_directories("results/spectrum");
          std::ofstream out("results/spectrum/gamma.txt");
          for (auto v : mergedHist_) out << v << "\n";
          G4cout << "[Run] Spectrum saved to results/spectrum/gamma.txt\n";
      }

      // Метод для злиття гістограм з потоків
      void MergeHistogram(const std::vector<int>& workerHist) {
          std::lock_guard<std::mutex> lock(histMutex_);
          if (mergedHist_.size() != workerHist.size()) {
              mergedHist_.resize(workerHist.size(), 0);
          }
          for (size_t i = 0; i < workerHist.size(); ++i) {
              mergedHist_[i] += workerHist[i];
          }
      }

      void RunCSV(G4String name) { /* як у тебе вже є — без змін */ }

      void ExportSpectrumTXT(G4String name) {
          if (!isMaster_) {
              G4cout << "ExportSpectrumTXT can only be called from the master thread." << G4endl;
              return;
          }
          std::filesystem::create_directories("results_txt");
          std::ofstream txt("results_txt/" + std::string(name) + "_spectrum.txt");
          for (size_t i = 0; i < mergedHist_.size(); ++i) {
              txt << mergedHist_[i] << "\n";
          }
          G4cout << "[Analysis] Simple spectrum written to results_txt/" << name << "_spectrum.txt\n";
      }
private:
    std::vector<int> mergedHist_;
    std::mutex histMutex_;

  private:
      G4GenericMessenger* messenger_ = nullptr;
      bool isMaster_ = false;
  };

  // ✅ Тільки AnalysisManager → ActionInitialization
  class ActionInitialization : public G4VUserActionInitialization {
  public:
      void Build() const override {
          if (!ENX09A::gParamsMessenger) ENX09A::gParamsMessenger = std::make_unique<ENX09A::ParamsMessenger>();
          SetUserAction(new PrimaryGeneratorAction());
          SetUserAction(new ProgressEventAction());
          SetUserAction(new AnalysisManager());
      }
      void BuildForMaster() const override {
          if (!ENX09A::gParamsMessenger) ENX09A::gParamsMessenger = std::make_unique<ENX09A::ParamsMessenger>();
          SetUserAction(new AnalysisManager());
      }
  };
} // namespace ENX09A
// ===== ПІДКЛЮЧЕННЯ ДО ЄДИНОГО МОСТУ ENTERPRISE =====
namespace Enterprise { int StartFromBridgeQt(int argc, char** argv); }
int main(int argc, char** argv) {
    G4cout << "===================================================" << G4endl;
    G4cout << " ENX09 Simulation Program Version 0.1" << G4endl;
    G4cout << " == Enterprise project Version 9.1 ==" << G4endl;
    G4cout << " Modeling particle interactions using a tungsten target" << G4endl;
    G4cout << " and a sensitive detector to record spectral data." << G4endl;
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
        }
    }
    runManager->SetUserInitialization(new ENX09A::MyDetector());
    runManager->SetUserInitialization(new ENX09A::MyPhysicsList());
    runManager->SetUserInitialization(new ENX09A::ActionInitialization());
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