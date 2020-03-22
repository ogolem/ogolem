pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                echo 'Building ogolem and running unit tests'
                sh 'gradle build -i'
            }
        }
        stage('Build ogolem manual') {
            steps {
                echo 'Building ogolem manual'
                sh 'cd manual && pdflatex manual.tex && bibtex manual && pdflatex manual.tex && pdflatex manual.tex && cd -'
            }
        }
        stage('Microbenchmark') {
            steps {
	        echo 'Running ogolem internal microbenchmarks'
                sh 'java -jar build/libs/ogolem-snapshot.jar --microbenchmarks'
            }
        }
        stage('Macrobenchmark - clusters') {
            steps {
	        echo 'Running ogolem internal macrobenchmarks for cluster structure optimization'
                sh 'java -jar build/libs/ogolem-snapshot.jar --macrobenchmarks -cluster cluster-macrobenchs.csv 2'
            }
        }
        stage('Macrobenchmark - adaptive') {
            steps {
	        echo 'Running ogolem internal macrobenchmarks for adaptive parameter optimization'
                sh 'java -jar build/libs/ogolem-snapshot.jar --macrobenchmarks -adaptive adaptive-macrobenchs.csv 2'
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
