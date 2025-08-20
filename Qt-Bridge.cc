#pragma once // Qt-Bridge-G4
#include <QApplication>
#include <G4UIExecutive.hh>
#include <G4VisExecutive.hh>
#include <G4UImanager.hh>

#include <thread>
#include <atomic>
#include <functional>
#include <iostream>

namespace Enterprise
{
    int StartFromBridge(int argc, char** argv) {
        std::string firstArg = (argc > 0) ? argv[0] : "none";
        G4cout << "Executable path: " << firstArg << G4endl;

        QApplication app(argc, argv);
        auto* ui = new G4UIExecutive(argc, argv);
        auto* visManager = new G4VisExecutive();
        visManager->Initialize();

        std::cout << "[QtBridge] Initializing Qt + Geant4 GUI\n";

        std::thread gui([&] {
            if (ui->IsGUI()) {
                ui->SessionStart();
                // ❌ НЕ встановлюємо running = false
            }
            });

        auto* uiMgr = G4UImanager::GetUIpointer();
        std::cout << "\n QtBridge interactive Geant4 interface\n"
            << " type Geant4 commands below, or 'exit' to quit\n";

        bool running = true;
        while (running) {
            std::cout << ">> ";
            std::string line;
            std::getline(std::cin, line);

            if (line == "exit") {
                running = false;
                break;
            }

            if (!line.empty()) {
                G4int code = uiMgr->ApplyCommand(line);
                if (code != 0)
                    std::cout << "[warn] Command failed with code " << code << "\n";
            }
        }

        gui.join();
        delete ui;
        delete visManager;
        return 0;
    }
}
