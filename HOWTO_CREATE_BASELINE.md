### How-to create a new spotbugs baseline.xml file

1) Comment
```
spotbugsMain {
  reports {
    html {
      enabled = true
      stylesheet = 'fancy-hist.xsl'
    }
  }
}
```
section in `build.gradle` to disable html output.

2) Remove old `baseline.xml`

3) Run `gradle.build` (it will fail)

4) `cp build/reports/spotbugs/main.xml baseline.xml`

5) Run `gradle.build` again (it should not fail now)
