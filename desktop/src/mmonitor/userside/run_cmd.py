from mmonitor.userside.MMonitorCMD import MMonitorCMD

def main():
    cmd_runner = MMonitorCMD()
    args = cmd_runner.parse_arguments()
    cmd_runner.initialize_from_args(args)
    cmd_runner.run()

if __name__ == "__main__":
    main()