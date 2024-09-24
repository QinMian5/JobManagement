# Author: Mian Qin
# Date Created: 9/23/24
import json


def main():
    X_STAR_list = [0, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1700, 1800]
    ramp_rate = 800 / 5000
    prd_time = 10000
    job_params = {}
    for X_STAR in X_STAR_list:
        ramp_time = int(X_STAR / ramp_rate)
        nsteps = int((ramp_time + prd_time) / 0.002)
        job_params[f"op_{X_STAR}"] = {
            "QBAR": {"TYPE": "parabola", "CENTER": X_STAR, "KAPPA": 0.05},
            "TEMPERATURE": 300,
            "RAMP_TIME": ramp_time,
            "PRD_TIME": prd_time,
            "NSTEPS": nsteps
        }
    with open("job_params.json", 'w') as file:
        json.dump(job_params, file, indent='\t')


if __name__ == "__main__":
    main()
