@echo off
REM PF 声学似然 A/B 验证入口（需 MATLAB R2020b+）
cd /d "%~dp0\.."
echo Working directory: %CD%
echo.
echo [1/3] 检查数据文件...
if not exist "Scenario\Step1_Data.mat" (
    echo 缺少 Step1_Data.mat，请先运行 Scenario\run_scenario.m
    pause
    exit /b 1
)
if not exist "Scenario\Step2_HeteroDetections.mat" (
    echo 缺少 Step2_HeteroDetections.mat，请先运行 step2_detection_simulation.m
    pause
    exit /b 1
)
echo.
echo [2/3] 运行 A/B 验证...
matlab -batch "run('validation/run_acoustic_likelihood_ab.m')"
if errorlevel 1 (
    echo MATLAB 执行失败。请确认 matlab 在 PATH 中。
    pause
    exit /b 1
)
echo.
echo [3/3] 完成。报告: report\PF_acoustic_likelihood_AB_验证报告.md
pause
