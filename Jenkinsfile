pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                echo 'Building ogolem and running unit tests'
                sh 'gradle build -i'
            }
        }
        stage('Microbenchmark') {
            steps {
	        echo 'Running ogolem internal microbenchmarks'
                sh 'java -jar build/libs/ogolem-snapshot.jar --microbenchmarks > microbenchmark_results.txt'
            }
        }
        stage('Build ogolem manual') {
            steps {
                echo 'Building ogolem manual'
                sh 'cd manual && pdflatex manual.tex && bibtex manual && pdflatex manual.tex && pdflatex manual.tex && cd -'
            }
        }
    }
    
    post {
        always {
            archiveArtifacts artifacts: 'build/libs/ogolem-snapshot.jar', fingerprint: true
            archiveArtifacts artifacts: 'manual/manual.pdf', fingerprint: true
            archiveArtifacts artifacts: 'microbenchmark_results.txt', fingerprint: true
        }
    }
}
