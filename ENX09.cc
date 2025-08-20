// ENX09М.cc — два послідовних Be-циліндра, джерело gamma з довільною енергією (макро або консоль),
// SD-диск на виході другого блоку, збір нейтронів (загальна кількість і спектр),
// багатопотокове злиття, прогрес у %, макро-керування, один файл.
//    Детектор:
//   - Корпус: поліетилен, прямокутний паралелепіпед 53×53×16,3 мм
//   - Чутливий об'єм (сцинтилятор): циліндр Ø39,30×8,55 мм, центр співпадає з центром куба, заповнений порошком
//   - Масова маса порошку: 31,64 мг
//   - Хімічна формула порошку: CoCO3*Co(OH)2*H2O (n=1)
//
#include "G4MTRunManager.hh"
#include "G4RunManager.hh"
#include "G4RunManagerFactory.hh"
#include "G4VUserDetectorConstruction.hh"
#include "G4VModularPhysicsList.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4EmExtraPhysics.hh"
#include "G4Neutron.hh"
#include "G4HadronPhysicsQGSP_BIC_HP.hh"
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
#include "G4ParticleGun.hh"
#include "G4Electron.hh"
#include "G4Gamma.hh"
#include "G4Positron.hh"
#include "G4Step.hh"
#include "G4HCofThisEvent.hh"
#include "G4PhysicsListHelper.hh"
#include "G4eMultipleScattering.hh"
#include "G4eIonisation.hh"
#include "G4eBremsstrahlung.hh"
#include "G4PhotoElectricEffect.hh"
#include "G4ComptonScattering.hh"
#include "G4GammaConversion.hh"
#include "G4Alpha.hh"
#include "G4Neutron.hh"
#include "G4Proton.hh"
#include "G4IonConstructor.hh"
#include "G4Decay.hh"
#include "G4RadioactiveDecay.hh"
#include "G4IonTable.hh"
#include "G4EmStandardPhysics.hh"
#include "G4DecayPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4HadronPhysicsFTFP_BERT.hh"
#include "G4VDiscreteProcess.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4DynamicParticle.hh"
#include "G4ParticleChange.hh"
#include "G4ProcessType.hh"
#include "G4Threading.hh"
#include "G4UImessenger.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4GenericMessenger.hh"
#include "Randomize.hh"
#include <mutex>
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
#include <locale>
#include <codecvt>
#include <cstdlib>
///  
namespace ENX09B
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
                ENX09B::BIN_WIDTH = binWidthCmd->GetNewDoubleValue(val) / CLHEP::MeV;
                ENX09B::Recompute();
            }
            else if (cmd == maxEnergyCmd) {
                ENX09B::MAX_EN = maxEnergyCmd->GetNewDoubleValue(val) / CLHEP::MeV;
                ENX09B::Recompute();
            }
            else if (cmd == gammaEnergyCmd) {
                ENX09B::PRIMARY_E = gammaEnergyCmd->GetNewDoubleValue(val) / CLHEP::MeV;
            }
        }
    private:
        G4UIcmdWithADoubleAndUnit* binWidthCmd{};
        G4UIcmdWithADoubleAndUnit* maxEnergyCmd{};
        G4UIcmdWithADoubleAndUnit* gammaEnergyCmd{};
        G4UIcmdWithABool* writeEntryCSVCmd{};
    };

    inline std::unique_ptr<ParamsMessenger> gParamsMessenger;
    // SD: ловити всі частинки у сцинтиляторі, веде два результуючі файли
    class ScintillatorSD : public G4VSensitiveDetector {
    public:
        ScintillatorSD(const G4String& name)
            : G4VSensitiveDetector(name) {
            std::filesystem::create_directories("results");
            spectrum_out_.open("results/scint_spectrum.csv", std::ios::out | std::ios::trunc);
            spectrum_out_ << "E_MeV;particle;isomeric;\n";
            cobalt_out_.open("results/cobalt_decay.csv", std::ios::out | std::ios::trunc);
            cobalt_out_ << "E_MeV;isomeric;\n";
        }
        ~ScintillatorSD() override {
            if (spectrum_out_.is_open()) spectrum_out_.close();
            if (cobalt_out_.is_open()) cobalt_out_.close();
        }
        G4bool ProcessHits(G4Step* step, G4TouchableHistory*) override {
            const auto* trk = step->GetTrack();
            const auto* def = trk->GetParticleDefinition();
            const auto* pre = step->GetPreStepPoint();
            double E = pre->GetKineticEnergy() / MeV;
            auto pos = pre->GetPosition();
            double x = pos.x() / mm, y = pos.y() / mm, z = pos.z() / mm;
            std::string pname = def->GetParticleName();
            bool isCobalt = (pname.find("Co") != std::string::npos);
            bool isomeric = false;
            // Перевірка на ізомерний стан (якщо ім'я містить 'm')
            if (isCobalt && pname.find("m") != std::string::npos) isomeric = true;
            // Запис у загальний спектр
            if (spectrum_out_.is_open()) {
                spectrum_out_ << ENX09B::fmt_comma(E) << ";" << pname << ";" << (isomeric ? 1 : 0) << ";"
                              << ENX09B::fmt_comma(x) << ";" << ENX09B::fmt_comma(y) << ";" << ENX09B::fmt_comma(z) << "\n";
            }
            // Запис у cobalt_decay.csv тільки для кобальту
            if (isCobalt && cobalt_out_.is_open()) {
                cobalt_out_ << ENX09B::fmt_comma(E) << ";" << (isomeric ? 1 : 0) << ";"
                            << ENX09B::fmt_comma(x) << ";" << ENX09B::fmt_comma(y) << ";" << ENX09B::fmt_comma(z) << "\n";
            }
            return true;
        }
    private:
        std::ofstream spectrum_out_;
        std::ofstream cobalt_out_;
    };

    class MyDetector : public G4VUserDetectorConstruction {
    public:
        G4VPhysicalVolume* Construct() override {
            auto* nist = G4NistManager::Instance();
            auto* vacuum = nist->FindOrBuildMaterial("G4_Galactic");
            auto* tungsten = nist->FindOrBuildMaterial("G4_W");
            auto* beryllium = nist->FindOrBuildMaterial("G4_Be");

            // --- ВСЯ ГЕОМЕТРІЯ ---
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
            // Маркер джерела (залишаємо для візуалізації)
            {
                auto* gunSolid = new G4Sphere("GunSolid", 0, 1.5 * mm, 0, 360 * deg, 0, 180 * deg);
                auto* gunLV = new G4LogicalVolume(gunSolid, vacuum, "GunLV");
                gunLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.7, 1.0)));
                new G4PVPlacement(nullptr, { 0,0,-5 * cm }, gunLV, "GunPV", worldLV, false, 0, true);
            }
            // Вольфрамова мішень (93×55×2 мм) у z = −5 см
            {
                auto* tgtSolid = new G4Box("Target", 93 * mm / 2, 55 * mm / 2, 1 * mm);
                auto* tgtLV = new G4LogicalVolume(tgtSolid, tungsten, "TargetLV");
                new G4PVPlacement(nullptr, { 0,0,-5 * cm }, tgtLV, "TargetPV", worldLV, false, 0, true);
                tgtLV->SetVisAttributes(wireOnly(G4Colour(0.4, 0.4, 0.4)));
            }
            // Перший Be-циліндр
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
            // SD-диск (залишаємо у геометрії, але НЕ робимо SD)
            const G4double SD_HALF_THICK = 0.25 * mm;
            const G4double sdZ = beZ2 + beThick / 2.0 + SD_HALF_THICK;
            auto* sdSolid = new G4Tubs("SDSolid", 0, beR, SD_HALF_THICK, 0, 360 * deg);
            auto* sdLV = new G4LogicalVolume(sdSolid, vacuum, "NeutronExitLV");
            sdLV->SetVisAttributes(wireOnly(G4Colour(1.0, 0.0, 0.0)));
            new G4PVPlacement(nullptr, { 0,0,sdZ }, sdLV, "NeutronExitPV", worldLV, false, 0, true);
            // === ДЕТЕКТОР === //
            const G4double scint_r = 39.30 / 2.0 * mm;
            const G4double scint_h = 8.55 * mm;
            const G4double detZ = beZ2 + beThick / 2.0 + 50.0 * mm;
            G4Element* elCo = nist->FindOrBuildElement("Co");
            G4Element* elC = nist->FindOrBuildElement("C");
            G4Element* elO = nist->FindOrBuildElement("O");
            G4Element* elH = nist->FindOrBuildElement("H");
            G4Material* scintMat = new G4Material("CoScint", 3.6 * g/cm3, 4);
            scintMat->AddElement(elCo, 2);
            scintMat->AddElement(elC, 1);
            scintMat->AddElement(elO, 5);
            scintMat->AddElement(elH, 2);
            auto* scintSolid = new G4Tubs("ScintSolid", 0, scint_r, scint_h/2, 0, 360*deg);
            fScintLV = new G4LogicalVolume(scintSolid, scintMat, "ScintLV");
            fScintLV->SetVisAttributes(wireOnly(G4Colour(0.2, 0.8, 0.2)));
            new G4PVPlacement(nullptr, {0,0,detZ}, fScintLV, "ScintPV", worldLV, false, 0, true);
            return worldPV;
        }
        void ConstructSDandField() override {
            auto* SDM = G4SDManager::GetSDMpointer();
            auto* sd = new ScintillatorSD("ScintillatorSD");
            SDM->AddNewDetector(sd);
            fScintLV->SetSensitiveDetector(sd);
        }
    private:
        G4LogicalVolume* fScintLV = nullptr;
    };
    // ── MyPhotoNuclearProcess
  // Цей клас реалізує процес фотоядерної взаємодії
    class MyPhotoNuclearProcess : public G4VDiscreteProcess {
    public:
        G4ParticleChange fParticleChange; // для вторинних частинок
        MyPhotoNuclearProcess() : G4VDiscreteProcess("MyPhotoNuclear") {
            SetProcessType(fHadronic);
            SetProcessSubType(131); // fPhotoNuclear
			SetVerboseLevel(1); // рівень виводу інформації
		}
        G4bool IsApplicable(const G4ParticleDefinition& particle) override {
            return (&particle == G4Gamma::Gamma());
        }

        G4VParticleChange* PostStepDoIt(const G4Track& aTrack, const G4Step& aStep) override {
            aParticleChange.Initialize(aTrack);

            const G4StepPoint* preStepPoint = aStep.GetPreStepPoint();
            const G4TouchableHandle& theTouchable = preStepPoint->GetTouchableHandle();
            G4VPhysicalVolume* thePhysical = theTouchable->GetVolume();
            if (!thePhysical) { return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep); }

            G4LogicalVolume* theLogical = thePhysical->GetLogicalVolume();
            if (theLogical->GetName() != "BeLV1" && theLogical->GetName() != "BeLV2") {
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
            }

            if (aTrack.GetKineticEnergy() < 1.67 * MeV) { // Threshold for 9Be(γ,n)8Be
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
            }

            double neutronEnergy = (aTrack.GetKineticEnergy() - 1.67 * MeV) * G4UniformRand();
            if (neutronEnergy <= 0) {
                return G4VDiscreteProcess::PostStepDoIt(aTrack, aStep);
            }

            aParticleChange.SetNumberOfSecondaries(1);
            G4DynamicParticle* neutron = new G4DynamicParticle(G4Neutron::Neutron(), aTrack.GetMomentumDirection(), neutronEnergy);
            aParticleChange.AddSecondary(neutron);
            aParticleChange.ProposeTrackStatus(fStopAndKill);

            return &aParticleChange;
        }

        // FIXED: Implemented the pure virtual function that was causing compile errors
        G4double GetMeanFreePath(const G4Track&, G4double, G4ForceCondition* cond) override
        {
            *cond = NotForced;
            return 5.0 * cm; // Fixed mean free path for simplicity
        }
    };
    // Фізичний список
    // (додаємо MyPhotoNuclearProcess для фотоядерних реакцій)
    class MyPhysicsList : public G4VModularPhysicsList {
    public:
        MyPhysicsList() {
            SetVerboseLevel(1);
            RegisterPhysics(new G4EmStandardPhysics());
            RegisterPhysics(new G4DecayPhysics());
            RegisterPhysics(new G4RadioactiveDecayPhysics());
            RegisterPhysics(new G4HadronPhysicsFTFP_BERT());
        }
        void ConstructProcess() override {
            G4VModularPhysicsList::ConstructProcess();
            // Додаємо кастомний PhotoNuclear процес для gamma
            auto* ph = G4PhysicsListHelper::GetPhysicsListHelper();
            ph->RegisterProcess(new MyPhotoNuclearProcess(), G4Gamma::Gamma());
        }
        void SetCuts() override {
            SetCutValue(0.001 * mm, "gamma");
            SetCutValue(0.001 * mm, "e-");
            SetCutValue(0.001 * mm, "e+");
            SetCutValue(0.001 * mm, "neutron");
            SetCutValue(0.001 * mm, "proton");
            SetCutValue(0.001 * mm, "alpha");
        }
    };
    class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction {
    public:
        PrimaryGeneratorAction() {
            electronGun_ = new G4ParticleGun(1);
            electronGun_->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("e-"));
            electronGun_->SetParticleEnergy(8.7 * MeV);
            electronGun_->SetParticleMomentumDirection({0,0,1});
            electronGun_->SetParticlePosition({0,0,-10 * cm});

            gammaGun_ = new G4ParticleGun(1);
            gammaGun_->SetParticleDefinition(G4ParticleTable::GetParticleTable()->FindParticle("gamma"));
            gammaGun_->SetParticleMomentumDirection({0,0,1});
            gammaGun_->SetParticlePosition({0,0,0}); // після вольфраму, перед берилієм
        }
        ~PrimaryGeneratorAction() override {
            delete electronGun_;
            delete gammaGun_;
        }

        void GeneratePrimaries(G4Event* evt) override {
            electronGun_->GeneratePrimaryVertex(evt);
            double energy = SampleGammaEnergy();
            gammaGun_->SetParticleEnergy(energy * MeV);
            gammaGun_->GeneratePrimaryVertex(evt);
        }
        // PDF для розподілу (Pade22)
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
        void SetEnergy(double en_MeV) { gammaGun_->SetParticleEnergy(en_MeV * MeV); }
    private:
        G4ParticleGun* electronGun_ = nullptr;
        G4ParticleGun* gammaGun_ = nullptr;
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
        }
        ~AnalysisManager() override { delete messenger_; }
        void EndOfRunAction(const G4Run*) override {
            if (isMaster_) {
                G4cout << "[Run] Дані зібрані ScintillatorSD (сцинтилятор).\n";
                G4cout << "[Run] Загальний спектр: results/scint_spectrum.csv\n";
                G4cout << "[Run] Розпади кобальта: results/cobalt_decay.csv\n";
            }
        }
    private:
        G4GenericMessenger* messenger_ = nullptr;
        bool isMaster_ = false;
    };

    // ✅ Тільки AnalysisManager → ActionInitialization
    class ActionInitialization : public G4VUserActionInitialization {
    public:
        void Build() const override {
            if (!ENX09B::gParamsMessenger) ENX09B::gParamsMessenger = std::make_unique<ENX09B::ParamsMessenger>();
            SetUserAction(new PrimaryGeneratorAction());
            SetUserAction(new ProgressEventAction());
            SetUserAction(new AnalysisManager());
        }
        void BuildForMaster() const override {
            if (!ENX09B::gParamsMessenger) ENX09B::gParamsMessenger = std::make_unique<ENX09B::ParamsMessenger>();
            SetUserAction(new AnalysisManager());
        }
    };
} // namespace ENX09B
// ===== ПІДКЛЮЧЕННЯ ДО ЄДИНОГО МОСТУ ENTERPRISE =====
namespace Enterprise { int StartFromBridgeQt(int argc, char** argv); }

int main(int argc, char** argv) {
    // Встановлюємо локаль для підтримки кирилиці у std::wcout/std::wcin
    std::locale loc("");
    std::locale::global(loc);
    std::wcin.imbue(loc);
    std::wcout.imbue(loc);

    G4cout << "===================================================" << G4endl;
    G4cout << " ENX09 Simulation Program Version 0.5" << G4endl;
    G4cout << " == Enterprise project Version 9.4 ==" << G4endl;
    G4cout << " Modeling gamma → (n,...) on two sequential Be blocks" << G4endl;
    G4cout << " Gamma source, variable energy, neutron spectrum after Be" << G4endl;
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
        else {
            G4cout << "[MT] dynamic_cast<G4MTRunManager*> failed! Not running in MT mode?\n";
        }
    }
    runManager->SetUserInitialization(new ENX09B::MyDetector());
    runManager->SetUserInitialization(new ENX09B::MyPhysicsList());
    runManager->SetUserInitialization(new ENX09B::ActionInitialization());
    // Ініціалізуєш сцену ДО запуску будь-якої макро або UI
    runManager->Initialize();
    if (threads == 0) {
        Enterprise::StartFromBridgeQt(argc, argv);
    }
    else {
        auto* UIm = G4UImanager::GetUIpointer();
        UIm->ApplyCommand("/control/execute run_phase.mac");
    }
    std::cout << "[INFO] Program finished. Press Enter to exit...\n";
    std::cin.get();
    delete runManager;
    return 0;
}