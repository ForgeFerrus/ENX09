cmake_minimum_required(VERSION 3.16)
project(ENX09 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_AUTOMOC ON)

# === Шлях до Qt: заміни на свою актуальну версію, якщо відрізняється
set(CMAKE_PREFIX_PATH "C:/Users/Forge/MyProject/Geant4/project_Qt/Qt-Bridge/lib/cmake")  # 🔁 замінити!

# === Qt модулі
find_package(Qt6 REQUIRED COMPONENTS Core Widgets)

# === Geant4
set(Geant4_DIR "C:/Users/Forge/MyProject/Geant4/program_files/lib/cmake/Geant4")  # 🔁 замінити!
find_package(Geant4 REQUIRED ui_all vis_all)

include_directories(${Geant4_INCLUDE_DIRS})
include_directories(
    "C:/Users/Forge/MyProject/Geant4/project_Qt/Qt-Bridge/include"
    "C:/Users/Forge/MyProject/Geant4/project_Qt/Qt-Bridge/include/QtWidgets"
)

# === Основна сцена Збірка цільового виконуваного файлу
add_executable(ENX09
    "C:/Users/Forge/MyProject/Geant4/project_Qt/Qt-Bridge/Qt-Bridge.cc"
    ENX09.cc
)

# === 🔗 Підключення бібліотек
target_link_libraries(ENX09
    PRIVATE Qt6::Core Qt6::Widgets
    ${Geant4_LIBRARIES}
)

# === 🛠 Додаткові include-и для Qt в релізі
target_include_directories(ENX09 PRIVATE
    ${Geant4_INCLUDE_DIRS}
    "C:/Users/Forge/MyProject/Geant4/project_Qt/Qt-Bridge/include"
    "C:/Users/Forge/MyProject/Geant4/project_Qt/Qt-Bridge/include/QtWidgets"
)