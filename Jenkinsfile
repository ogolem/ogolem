pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                echo 'Building  ogolem and running unit tests'
                sh 'gradle build -i'
            }
        }
    }

    post {
        always {
            archiveArtifacts artifacts: 'build/libs/ogolem-snapshot.jar', fingerprint: true
            junit 'build/reports/tests/test/*'
        }
    }
}
