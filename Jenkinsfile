pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                echo 'Building ogolem and running unit tests'
                sh 'gradle build -i -Dorg.gradle.java.home=/usr/local/openjdk16'
                // touch test reports to avoid junit complaining about old, cached ones
                sh 'find . -name "TEST-*.xml" -exec touch {} \\;'
                junit 'build/test-results/**/*.xml'
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
                sh '/usr/local/openjdk16/bin/java --add-modules jdk.incubator.vector -jar build/libs/ogolem-snapshot.jar --microbenchmarks'
            }
        }
        stage('Macrobenchmark - clusters') {
            steps {
	        echo 'Running ogolem internal macrobenchmarks for cluster structure optimization'
                sh '/usr/local/openjdk16/bin/java --add-modules jdk.incubator.vector -jar build/libs/ogolem-snapshot.jar --macrobenchmarks -cluster cluster-macrobenchs.csv 2'
            }
        }
        stage('Macrobenchmark - adaptive') {
            steps {
	        echo 'Running ogolem internal macrobenchmarks for adaptive parameter optimization'
                sh '/usr/local/openjdk16/bin/java --add-modules jdk.incubator.vector -jar build/libs/ogolem-snapshot.jar --macrobenchmarks -adaptive adaptive-macrobenchs.csv 2'
            }
        }
    }
    
    post {
        always {
            archiveArtifacts artifacts: 'build/libs/ogolem-snapshot.jar', fingerprint: true
            archiveArtifacts artifacts: 'manual/manual.pdf', fingerprint: true
            archiveArtifacts artifacts: 'build/reports/spotbugs/main.html', fingerprint: true
        }
    }
}
