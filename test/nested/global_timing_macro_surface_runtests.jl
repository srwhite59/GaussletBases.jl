@testset "Global timing macro surface" begin
    old_config = GaussletBases.TimeG._TIMING_CONFIG[]
    try
        @test timing_enabled() == GaussletBases.TimeG.timing_enabled()
        @test timing_live_enabled() == GaussletBases.TimeG.timing_live_enabled()

        reset_timing_report!()
        set_timing!(true)
        set_timing_live!(false)
        set_timing_thresholds!(expand = 0.0, drop = 0.0)

        @timeg "outer" begin
            sleep(0.002)
            @timeg "inner" begin
                sleep(0.001)
            end
        end

        report = current_timing_report()
        @test report isa GaussletBases.TimeG.TimingReport
        @test length(report.roots) == 1
        root = only(report.roots)
        @test root.label == "outer"
        @test root.elapsed_seconds > 0.0
        @test root.self_seconds >= 0.0
        @test root.call_count == 1
        @test length(root.children) == 1
        child = only(root.children)
        @test child.label == "inner"
        @test child.elapsed_seconds > 0.0

        rendered = timing_report(report)
        @test occursin("GaussletBases timing report", rendered)
        @test occursin("outer", rendered)
        @test occursin("inner", rendered)

        live_path, live_io = mktemp()
        close(live_io)
        try
            open(live_path, "w") do io
                redirect_stdout(io) do
                    reset_timing_report!()
                    set_timing!(true)
                    set_timing_live!(true)
                    set_timing_thresholds!(expand = 0.0, drop = 0.0)
                    @timeg "live outer" begin
                        @timeg "live inner" begin
                            sleep(0.001)
                        end
                    end
                end
            end
            live_output = read(live_path, String)
            @test occursin("live outer: ", live_output)
            @test occursin("live inner: ", live_output)
            @test occursin("seconds", live_output)
        finally
            rm(live_path; force = true)
        end

        reset_timing_report!()
        set_timing!(false)
        set_timing_live!(false)
        @timeg "disabled" begin
            sleep(0.001)
        end
        disabled_report = current_timing_report()
        @test isempty(disabled_report.roots)
    finally
        GaussletBases.TimeG._TIMING_CONFIG[] = old_config
        reset_timing_report!()
    end
end
