
// ENX09C.cc — фізична модель детектора згідно експериментальної задачі
//
// Детектор:
//   - Корпус: поліетилен, прямокутний паралелепіпед 53×53×16,3 мм
//   - Порожнина для порошку: циліндр Ø39,30×8,55 мм, центр співпадає з центром куба
//   - Масова маса порошку: 31,64 мг
//   - Хімічна формула порошку: CoCO3*Co(OH)2*nH2O
//
// Джерело:
//   - Енергія лінії 60mCo: 58,59 кеВ
//   - Енергія ліній 60Co: 1173,2 кеВ, 1332,5 кеВ
//
// Додатково (для аналізу):
//   - Енергії іонізації (1с, 1s2) для детекторних матеріалів:
//     32Ge: 14,119 кеВ / 13,557 кеВ
//     11Na: 1,649 кеВ / 1,465 кеВ
//     53I: 39,721 кеВ / 38,717 кеВ
//     81Tl: 98,592 кеВ / 96,783 кеВ
//
// Мета: моделювання проходження гамма-квантів та вторинних частинок через фізичну мішень, збір спектру та просторового розподілу у відповідності до реальної конструкції.
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
#include "G4hMultipleScattering.hh"
#include "G4hIonisation.hh"
#include "G4HadronElasticProcess.hh"
#include "G4NeutronCaptureProcess.hh"
#include "G4NeutronFissionProcess.hh"
#include "G4ionIonisation.hh"
#include "G4ProcessManager.hh"
#include "G4CascadeInterface.hh"
#include "G4LeptonConstructor.hh"
// Моделі та XS-таблиці
#include "G4HadronElastic.hh"
#include "G4NeutronElasticXS.hh"
#include "G4NeutronRadCapture.hh"
#include "G4NeutronCaptureXS.hh"
#include "G4LFission.hh"
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

///  
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
    // SD: ловить gamma перед Be; веде гістограму і лог хітів у CSV per-thread
    class GammaEntrySD : public G4VSensitiveDetector {
    public:
        // Додаємо метод для доступу до гістограми (має бути public)
        const std::vector<int>& GetHistogram() const { return hist_; }
        GammaEntrySD(const G4String& name)
            : G4VSensitiveDetector(name), hist_(ENX09C::NUM_BINS, 0) {
            // Завжди створюємо results і відкриваємо CSV
            std::filesystem::create_directories("results");
            static std::mutex mtx;
            std::lock_guard<std::mutex> lock(mtx);
            static bool opened = false;
            if (!opened) {
                out_.open("results/neutron_entry.csv", std::ios::out | std::ios::trunc);
                out_ << "E_MeV;x_mm;y_mm;z_mm;ux;uy;uz;theta_deg;r_mm\n";
                opened = true;
            } else {
                out_.open("results/neutron_entry.csv", std::ios::out | std::ios::app);
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
            if (trk->GetParticleDefinition()->GetParticleName() != "neutron") return false;
            const auto* pre = step->GetPreStepPoint();
            const double E = pre->GetKineticEnergy() / MeV;
            const int nb = static_cast<int>(hist_.size());
            int bin = static_cast<int>(E / ENX09C::BIN_WIDTH);
            if (bin >= 0 && bin < nb) hist_[bin]++;

            if (out_.is_open()) {
                const auto pos = pre->GetPosition();
                const double x = pos.x() / mm, y = pos.y() / mm, z = pos.z() / mm;
                const auto dir = pre->GetMomentumDirection();
                const double ux = dir.x(), uy = dir.y(), uz = dir.z();
                static constexpr double pi = 3.14159265358979323846;
                const double theta = std::acos(std::clamp(uz, -1.0, 1.0)) * 180.0 / pi;
                const double r = std::hypot(x, y);
                out_ << ENX09C::fmt_comma(E) << ";"
                    << ENX09C::fmt_comma(x) << ";"
                    << ENX09C::fmt_comma(y) << ";"
                    << ENX09C::fmt_comma(z) << ";"
                    << ENX09C::fmt_comma(ux) << ";"
                    << ENX09C::fmt_comma(uy) << ";"
                    << ENX09C::fmt_comma(uz) << ";"
                    << ENX09C::fmt_comma(theta) << ";"
                    << ENX09C::fmt_comma(r) << "\n";
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
            auto* beryllium = nist->FindOrBuildMaterial("G4_Be");
            auto* polyethylene = nist->FindOrBuildMaterial("G4_POLYETHYLENE");

            // Додаємо матеріал порошку CoCO3*Co(OH)2*nH2O (приблизна формула)
            // Молярна маса ~ 209.9 г/моль (CoCO3) + 92.94 г/моль (Co(OH)2) + n*18 г/моль (H2O)
            // Для простоти: беремо CoCO3 + Co(OH)2 + 2H2O = 320.84 г/моль
            // Густина: маса/об'єм = 31.64 мг / (π*(19.65 мм)^2*8.55 мм) [мм^3]
            const G4double powder_mass_mg = 31.64; // мг
            const G4double powder_mass_g = powder_mass_mg * 1e-3; // г
            const G4double holeR = 19.65 * mm; // O39.3 мм
            const G4double holeH = 8.55 * mm / 2.0; // половина висоти
            const G4double powder_vol_mm3 = M_PI * holeR * holeR * (2*holeH); // мм^3
            const G4double powder_vol_cm3 = powder_vol_mm3 * 1e-3; // см^3
            const G4double powder_density = powder_mass_g / powder_vol_cm3; // г/см^3

            G4Element* elCo = nist->FindOrBuildElement("Co");
            G4Element* elC = nist->FindOrBuildElement("C");
            G4Element* elO = nist->FindOrBuildElement("O");
            G4Element* elH = nist->FindOrBuildElement("H");
            G4Material* powder = new G4Material("CoPowder", powder_density * g/cm3, 4);
            powder->AddElement(elCo, 2);
            powder->AddElement(elC, 1);
            powder->AddElement(elO, 7); // CO3 + (OH)2 + 2H2O = 3+2+2=7 O
            powder->AddElement(elH, 4); // (OH)2 + 2H2O = 2+2=4 H

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
            auto solidRed = [](G4Colour c) {
                auto* vis = new G4VisAttributes(c);
                vis->SetForceSolid(true);
                vis->SetForceWireframe(false);
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

            // SD-кубічний детектор (корпус поліетилен)
            const G4double detX = 53.0 * mm;
            const G4double detY = 53.0 * mm;
            const G4double detZ = 16.3 * mm;
            const G4double detCenterZ = 316.4 * mm; // 31.64 см
            auto* detSolid = new G4Box("CubeDet", detX/2, detY/2, detZ/2);
            fGammaEntryLV = new G4LogicalVolume(detSolid, polyethylene, "CubeDetLV");
            fGammaEntryLV->SetVisAttributes(solidRed(G4Colour(1.0, 0.0, 0.0)));
            auto* detPV = new G4PVPlacement(nullptr, { 0,0,detCenterZ }, fGammaEntryLV, "CubeDetPV", worldLV, false, 0, true);

            // Порожнина-мішень (циліндр, заповнений порошком)
            auto* holeSolid = new G4Tubs("CubeHole", 0, holeR, holeH, 0, 360 * deg);
            auto* holeLV = new G4LogicalVolume(holeSolid, powder, "CubeHoleLV");
            holeLV->SetVisAttributes(wireOnly(G4Colour(1.0, 1.0, 1.0)));
            // Центр отвору співпадає з центром куба
            new G4PVPlacement(nullptr, {0,0,0}, holeLV, "CubeHolePV", fGammaEntryLV, false, 0, true);

            // Додатковий детектор (ще один куб після всього, на 3 см далі по z)
            const G4double det2CenterZ = detCenterZ + detZ/2 + 30.0 * mm + detZ/2; // 3 см після першого
            auto* det2Solid = new G4Box("CubeDet2", detX/2, detY/2, detZ/2);
            auto* det2LV = new G4LogicalVolume(det2Solid, polyethylene, "CubeDet2LV");
            det2LV->SetVisAttributes(solidRed(G4Colour(0.0, 0.7, 0.0))); // зелений solid
            new G4PVPlacement(nullptr, {0,0,det2CenterZ}, det2LV, "CubeDet2PV", worldLV, false, 0, true);

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
            if (track.GetKineticEnergy() < 7.0 * MeV)
            return &fParticleChange;

        // Створюємо нейтрон як простий вторинний
        auto* neutron = new G4DynamicParticle(G4Neutron::Neutron(), RandomDirection(), 5.0 * MeV);
        fParticleChange.SetNumberOfSecondaries(1);
        fParticleChange.AddSecondary(neutron);

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
    MyPhysicsList() {
        SetVerboseLevel(1); // рівень виводу інформації
        G4IonConstructor().ConstructParticle(); // конструктор іонів
    }
    void ConstructParticle() override {
        // Створення частинок
        G4Electron::ElectronDefinition();
        G4Gamma::GammaDefinition();
        G4Positron::PositronDefinition();  // ← ось вона — ключ до тиші!
		G4Proton::Proton();
		// Створення нейтронів та альфа-частинок
        G4Neutron::Neutron();
        G4Alpha::Alpha();
		// Створення всіх іонів

    }
    void ConstructProcess() override {
        AddTransportation(); // 🔑 обов’язково перший

        auto* helper = G4PhysicsListHelper::GetPhysicsListHelper();

        // електрони
        auto* eMinus = G4Electron::Electron();
        helper->RegisterProcess(new G4eMultipleScattering(), eMinus);
        helper->RegisterProcess(new G4eIonisation(), eMinus);
        helper->RegisterProcess(new G4eBremsstrahlung(), eMinus);

        // позитрони
        auto* ePlus = G4Positron::Positron();
        helper->RegisterProcess(new G4eMultipleScattering(), ePlus);
        helper->RegisterProcess(new G4eIonisation(), ePlus);
        helper->RegisterProcess(new G4eBremsstrahlung(), ePlus);

        // гамма
        auto* gamma = G4Gamma::Gamma();
        helper->RegisterProcess(new G4PhotoElectricEffect(), gamma);
        helper->RegisterProcess(new G4ComptonScattering(), gamma);
        helper->RegisterProcess(new G4GammaConversion(), gamma);
        // наш кастомний γ → n процес
        auto* myPhoto = new MyPhotoNuclearProcess(); // з класу вище
        helper->RegisterProcess(myPhoto, gamma);

        // альфа
        auto* alpha = G4Alpha::Alpha();
        helper->RegisterProcess(new G4hMultipleScattering(), alpha);
        helper->RegisterProcess(new G4ionIonisation(), alpha);

        // протон
        auto* p = G4Proton::Proton();
        helper->RegisterProcess(new G4hMultipleScattering(), p);
        helper->RegisterProcess(new G4ionIonisation(), p);

        // нейтрон
        auto* n = G4Neutron::Neutron();
        auto* elastic = new G4HadronElasticProcess();
        elastic->RegisterMe(new G4HadronElastic());
        elastic->AddDataSet(new G4NeutronElasticXS());
        helper->RegisterProcess(elastic, n);

        auto* capture = new G4NeutronCaptureProcess();
        capture->RegisterMe(new G4NeutronRadCapture());
        capture->AddDataSet(new G4NeutronCaptureXS());
        helper->RegisterProcess(capture, n);

        auto* fission = new G4NeutronFissionProcess();
        fission->RegisterMe(new G4LFission());
        helper->RegisterProcess(fission, n);
    }

    void SetCuts() override {
        for (const auto& name : { "gamma", "e-", "e+", "neutron", "proton", "alpha" })
            SetCutValue(0.001 * mm, name);
    }
};

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
    // ✅ AnalysisManager — єдиний RunAction
    class AnalysisManager : public G4UserRunAction {
  public:
      AnalysisManager() {
          isMaster_ = !G4Threading::IsWorkerThread();
          messenger_ = new G4GenericMessenger(this, "/analysis/", "Analysis control");
          if (isMaster_) {
              messenger_->DeclareMethod("run", &AnalysisManager::RunCSV)
                  .SetGuidance("Merge spectrum + generate CSV").SetParameterName("name", false);
              // Видалено exportSpectrumTXT
          }
      }
      ~AnalysisManager() override { delete messenger_; }

      void EndOfRunAction(const G4Run*) override {
          if (!isMaster_) {
              auto* base = G4SDManager::GetSDMpointer()->FindSensitiveDetector("GammaEntrySD");
              if (auto* sd = dynamic_cast<ENX09C::GammaEntrySD*>(base)) {
                  const auto& localHist = sd->GetHistogram();
                  auto* master = const_cast<AnalysisManager*>(static_cast<const AnalysisManager*>(G4MTRunManager::GetMasterRunManager()->GetUserRunAction()));
                  if (master) master->MergeHistogram(localHist);
              }
              return;
          }
          // Зберігаємо тільки CSV спектр
          std::filesystem::create_directories("results/spectrum");
          std::ofstream csv("results/spectrum/neutron_spectrum.csv");
          csv << "E_MeV;count\n";
          for (size_t i = 0; i < mergedHist_.size(); ++i) {
              double e = (i + 0.5) * ENX09C::BIN_WIDTH;
              csv << ENX09C::fmt_comma(e) << ";" << mergedHist_[i] << "\n";
          }
          G4cout << "[Run] Neutron spectrum saved to results/spectrum/neutron_spectrum.csv\n";
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

    void RunCSV(G4String name) { /* залишено для сумісності, нічого не робить */ }
    // Видалено всі методи, що експортують TXT-файли
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
    G4cout << " ENX09C Simulation Program Version 1.0" << G4endl;
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